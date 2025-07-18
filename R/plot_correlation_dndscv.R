#' Plot the correlation of synonymous mutations to observed mutations
#'
#' @description
#' Calculate and plot the correlation of the synonymous mutations to the observed mutations
#' as a quality measure of the gene-specific mutation model to assess the addition of different
#' measures to improve the correlation. Base is the length of the gene, additions are either the trinucleotides
#' and the covariates. For correlation metric, pearson correlation is used.
#' #'
#' @param dnds dnds output list resulting from \link[dndscv::dndscv]{dndscv}
#' @param L 3-dimensional array for trinucleotides per gene
#'
#' @returns Plot with the correlation values
#' @export
#'
#' @examples
plot_correlation_dnsdcv = function(dnds, L) {

  avg_rate = sum(dnds$N)/sum(dnds$L)
  gm = dnds$genemuts # shorthand for genemuts
  ratio_cv = gm$exp_syn_cv / gm$exp_syn

  df_syn = data.frame(
    observed_synonymous = gm$n_syn,
    length = colSums(L[,1,]) * avg_rate,
    `length + trinuc` = gm$exp_syn,
    `length + trinuc + covs` = gm$exp_syn_cv)

    data.frame(correlation = cor(df_syn)[-1,1]) |>
      mutate(name = colnames(df_syn)[-1]) |>
      ggplot(aes(x = name, y = correlation)) +
      geom_col() +
      theme_classic() +
      scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
      geom_text(aes(label = format(round(correlation, 2), trim = T, drop0trailing = TRUE, big.mark = ","),
                    y = correlation), vjust = -0.2, position = position_dodge(0.9)) +
      labs(x = NULL)
}

#' Plot the correlation of synonymous mutations to observed mutations
#'
#' @description
#' Calculate and plot the correlation of the synonymous mutations to the observed mutations
#' as a quality measure of the gene-specific mutation model to assess the addition of different
#' measures to improve the correlation. Base is the lenght of the gene, additions are either the trinucleotides
#' and the covariates. For correlation metric, pearson correlation is used.
#' #'
#' @param dnds_intron Object resulting from \link[wintr::dndscv_intron]{dndscv_intron}
#' @param L 3-dimensional array for trinucleotides per gene
#'
#' @returns Plot with the correlation values
#' @export
#'
#' @examples
plot_correlation_dnsdcv_intron = function(dnds_intron, L) {

  avg_rate = sum(dnds_intron$N)/sum(dnds_intron$L)
  gm = dnds_intron$genemuts # shorthand for genemuts
  ratio_cv = gm$exp_syn_cv / gm$exp_syn

  wintr = dnds_intron$globaldnds[6,2]

  length = (colSums(L[,1,]) + (colSums(L[,5,]) *wintr))
  df_syn = data.frame(
    observed_synonymous = gm$n_syn + gm$n_intr,
    length = length * avg_rate,
    `length + covs` = length * avg_rate * ratio_cv,
    `length + trinuc` = gm$exp_syn + gm$exp_intr,
    `length + trinuc + covs` = gm$exp_intr_syn_cv)

  data.frame(correlation = cor(df_syn)[-1,1]) |>
    mutate(name = colnames(df_syn)[-1]) |>
    ggplot(aes(x = name, y = correlation)) +
    geom_col() +
    theme_classic() +
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
    geom_text(aes(label = format(round(correlation, 2), trim = T, drop0trailing = TRUE, big.mark = ","),
                  y = correlation), vjust = -0.2, position = position_dodge(0.9)) +
    labs(x = NULL)
}
