#' plot correlation rate for 1, 12, 96 and 192 trinucleotide rates across the individual samples
#'
#' @param dnds_intron dnds output list resulting from \link[dndscv::dndscv]{dndscv} or
#'  \link[wintr::dndscv_intron]{dndscv_intron}
#'
#' @returns Plot with the correlation values for
#' @export
#'
#' @examples
plot_correlation_substmodel = function(dnds_intron) {

  # get the average mutation rate compared against trinucs:
  df_syn_rates = triplet_match_substmodel
  flat_rate = sum(dnds_intron$N)/sum(dnds_intron$L)

  # observed_rates- (mean of intron + syn)

  N = rowSums(dnds_intron$N, dims = 2)
  L = rowSums(dnds_intron$L, dims = 2)
  muts_syn = N[,1]
  muts_intr = N[,5] / dnds_intron$globaldnds[6,2]
  obs_rates = (muts_syn + muts_intr) / L[,5]

  # 192 submodel rates
  rates_192 = get_syn_mutrates(dnds_intron)

  # 12 base mutation rates
  submodel = dnds_intron$mle_submodel_12
  rates =  setNames(submodel[,2], submodel[,1])
  rates_12 = sapply(substmodel_introns_12[,1], function(x) prod(rates[base::strsplit(x,split="\\*")[[1]]]))
  # Expected rate per available site

  rate_model = df_syn_rates |>
    mutate(obs_rates = obs_rates,
           `flat (1)` = flat_rate,
           `trinucs-tx (192)` = rates_192,
           `subs-tx (12)` = rates_12) |>
    group_by(triplet) |>  mutate(`trinucs (96)` = mean(`trinucs-tx (192)`)) |>
    group_by(type) |>  mutate(`subs (6)` = mean(`subs-tx (12)`))  |>
    ungroup() |>
    select(-strand, -trinuc, -mut, triplet, -type, -triplet)  |>
    column_to_rownames("mut_type")

  rate_model_diff = rate_model |>
    mutate(across(everything(),  ~ . / obs_rates)) |>
    rownames_to_column("mut_type") |>
    left_join(triplet_match_substmodel) |>
    dplyr::select(-obs_rates)

  model_difference = rate_model_diff |>
    pivot_longer(-c(mut_type, strand, trinuc, mut, triplet, type)) |>
    mutate(name = as.factor(name))

    ggplot(model_difference, aes(x = name, y = value, color = name)) +
      geom_boxplot() +
      ggbeeswarm::geom_quasirandom() +
      geom_hline(yintercept = 1, linetype = "dotted") +
      facet_wrap(. ~ type) +
      theme_classic() +
      cowplot::panel_border() +
      ggsci::scale_color_igv() +
      labs(y = "Relative rate compared with\naverage observed 192-trinucleotide rate",
           x = "mutation model (number of categories)") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
}

