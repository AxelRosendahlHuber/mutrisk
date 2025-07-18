library(patchwork)

#' Get Synonymous Mutation Rates
#'
#' This function calculates the expected synonymous mutation rates per available
#' site based on a provided dN/dS submodel.
#'
#' @param dnds A list containing dN/dS analysis results, specifically expected to have a
#'   `mle_submodel` component. `mle_submodel` should be a data frame or matrix
#'   where the first column (`[,1]`) contains mutation types (e.g., "A*C", "T*G")
#'   and the second column (`[,2]`) contains their corresponding rates.
#' @return A named numeric vector where names are the mutation types (e.g., "A*C")
#'   and values are their calculated expected rates per available site.
get_syn_mutrates = function(dnds) {
  submodel = dnds$mle_submodel
  rates =  setNames(submodel[,2], submodel[,1])
  sapply(substmodel[,1], function(x) prod(rates[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
}

#' Compare Intron and Exon dN/dS and Mutation Rates
#'
#' This function performs a comprehensive comparison of mutation patterns and
#' selection parameters between exon and intron regions based on dN/dS analysis
#' results. It generates a multi-panel plot summarizing the differences.
#'
#' @param dnds_exon A list containing dN/dS analysis results for exon regions.
#'   Expected to have `nsites` (data frame with `n_muts` and `n_sites`) and
#'   `globaldnds` (data frame with `name`, `mle`, `cilow`, `cihigh`) components,
#'   and a structure compatible with `get_syn_mutrates`.
#' @param dnds_intron A list containing dN/dS analysis results for intron regions.
#'   Expected to have `nsites` (data frame with `n_muts` and `n_sites`) and
#'   `globaldnds` (data frame with `name`, `mle`, `cilow`, `cihigh`) components,
#'   and a structure compatible with `get_syn_mutrates`.
#' @return A `patchwork` object representing a multi-panel plot, comparing:
#'   1. Number of mutations and mutable sites between exon and intron-exon regions.
#'   2. Global dN/dS ($\Omega$) estimates and their confidence intervals.
#'   3. Trinucleotide mutation rates (absolute and exon/intron ratio).
#'   4. A log-log scatter plot of exon vs. intron mutation rates.
#' @importFrom dplyr mutate left_join filter
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap theme scale_y_continuous geom_text labs geom_errorbar geom_hline theme_classic scale_fill_manual scale_alpha_manual geom_point geom_abline scale_y_log10 scale_x_log10
#' @importFrom cowplot theme_cowplot
#' @importFrom data.table data.table rbindlist
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#'
#' @export
compare_intron_exon = function(dnds_exon, dnds_intron) {

  counts = data.frame(
    type = c("number of mutations", "number of mutable sites"),
    exon= c(sum(dnds_exon$nsites$n_muts),
                    sum(dnds_exon$nsites$n_sites) / 3),
    intron_exon = c(sum(dnds_intron$nsites$n_muts),
              sum(dnds_intron$nsites$n_sites)/ 3 )) |>
    mutate(exon = as.integer(exon),
           intron_exon = as.integer(intron_exon),
      ratio = round(intron_exon/exon, 1)) |>
    mutate(plotname = paste0(type, "\nRatio: 1/", ratio))

  counts_plot = counts |>
    pivot_longer(-c(type, ratio, plotname)) |>
    ggplot(aes(x = name, y = value, fill = name)) +
    facet_wrap(.~ plotname  , scales = "free_y") +
    geom_col() +
    cowplot::theme_cowplot() +
    theme(legend.position = "none") +
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
    geom_text(aes(label = format(value, trim = T, drop0trailing = TRUE, big.mark = ","), y = value), vjust = -0.2, position = position_dodge(0.9)) +
    labs(x = NULL, y = "number of sites/mutations")

  # compare differences in selection parameters
  omega_plot = rbindlist(list(dnds_exon = dnds_exon$globaldnds, dnds_intron = dnds_intron$globaldnds), idcol = "model") |>
    ggplot(aes(x = name, fill = model, y = mle, ymin = cilow, ymax = cihigh)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(position = position_dodge(width = 0.9), width = 0.5) +
    geom_hline(yintercept = 1) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult=c(0,0.1)), limits = c(0, 2.5)) +
    labs(x = NULL, y = expression(Omega ~ "  estimated")) +
    theme(legend.position = "none")

  mutrates_exon = get_syn_mutrates(dnds_exon)
  mutrates_intron = get_syn_mutrates(dnds_intron)

  df = data.table(mut_type = names(mutrates_exon), `mutrates exon` = mutrates_exon,
                  `mutrates intron` = mutrates_intron) |>
    left_join(triplet_match_substmodel, by = "mut_type") |>
    mutate(`exon/intron ratio` = mutrates_exon / mutrates_intron)

  df_long = df |>
    pivot_longer(c(`mutrates exon`, `mutrates intron`, `exon/intron ratio`))

  trinuc_plot = ggplot(df_long, aes(x = triplet, y = value, alpha = strand, fill = type)) +
    geom_col(position = position_dodge()) + facet_grid(name ~ . , scales = "free_y") +
    geom_hline(data = filter(df_long, name == 'exon/intron ratio'), aes(yintercept = 1), linetype = "dashed", color = "black") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = COLORS6) +
    scale_alpha_manual(values = c(0.5, 1)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = NULL, y = "(relative) mutrate", fill = NULL)

  df_log_mutrates = df |>
    ggplot(aes(x = `mutrates exon`, y = `mutrates intron`, color = type)) +
    geom_point() +
    geom_abline(slope = 1, linetype = "dotted") +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(values = COLORS6) +
    theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.7),
          legend.background = element_blank()) +
    labs(color = NULL)


  layout = "
  AADD
  BCDD"
  mg = 5
  counts_plot + omega_plot + df_log_mutrates + trinuc_plot +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = "A") &
    theme(plot.margin = margin(mg, mg, mg, mg, unit = "mm"))
}
