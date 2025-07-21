# Functions to calculate the VAF and sensitivity

# Estimate the 'true' underlying mean variant allele fraction (VAF) from a list of VAFs from a set of mutations
#' @description
#' Performs a density using the \link[stats]{stats} package \link[stats::density]{density} function.
#' The estimated 'true' VAF is the estimated maximum of the density curve. When estimated VAF >0.5 it will be set to 0.5
#'
#' @param vaf_values Numeric vector of the observed VAF of the mutations
#'
#' @returns Numeric value estimating the VAF
#' @export
#'
#' @examples
estimate_vaf = function(vaf_values) {
  density_est <- density(vaf_values, bw = 0.1) # draw a density curve for the values
  vaf = density_est$x[which.max(density_est$y)] # find the VAF as the max of the density curve
  if (vaf > 0.5) {
    vaf = 0.5
  } # biologically, the VAF cannot be higher than 0.5
  return(vaf)
}

#' Calculate the sensitivity to detect a mutation, assuming the occurrence of mutations is Poisson distributed.
#'
#' @description
#' Function to calculate the sensitivity to detect a mutation, using the
#' estimated mean variant allele fraction (vaf) (can be estimated using \code{\link{estimate_vaf}}), the depth.
#' In addition, the `minalt` parameter can be used to filter out samples with too low
#'
#' @param coverage Numeric vector of the mean coverage of the sample
#' @param vaf Numeric vector with the estimated variant allele fraction (VAF)
#' @param minalt minimum alternative allele count needed to register a mutation
#'
#' @returns a vector of the estimated sensitivity, ranging from 0 (zero sensitivity) to 1 (max sensitivity)
#' @export
#'
#' @examples
get_sensitivity <- function(coverage, vaf, minalt = 4) {
  mapply(
    function(cov, v) {
      sum(rpois(1e5, lambda = cov * v) >= minalt) / 1e5
    },
    coverage,
    vaf
  )
}

#' Test for the effect of coverage on the number of observed mutations
#'
#' Function performs regression analysis using a linear mixed model of the mutations to the coverage.
#' Likelihood ratio tests are performed to
#'
#' @param cell_muts data frame with mutations, containing
#' @param metadata
#'
#' @returns
#' a list with two elements:
#' 1. df_results: dataframe with the p-values for the likelihood ratio tests of the lmmm's testing for the effect of coverage.
#' - unadj: unadjusted mutations are used
#' - adj: mutations are corrected (divided by) the sensitivity
#' 2 plot indicating the raw values
#' @export
#'
#' @examples
effect_coverage_vaf = function(cell_muts, metadata) {
  # initialize empty dataframe and plot list
  df_results = data.frame(pval_coverage_unadj = NA,
                          pval_coverage_adj = NA)

  plot_list = list()

  # peform the analysis group by group
  muts_count = cell_muts[ , .N, by = "sampleID"]
  df_mut_count = inner_join(muts_count, metadata, by = "sampleID") |>
      dplyr::mutate(coverage = as.numeric(coverage),
                    Nadj = N / sensitivity)

  lmm_model = lme4::lmer(N ~ age + coverage + (1 | donor),
                         data = df_mut_count)
  lmm_null <- lme4::lmer(N ~ age + (1 | donor), data = df_mut_count) # Null model without predictors
  anv = anova(lmm_null, lmm_model) # Likelihood ratio test
  df_results[ , 1] = anv$`Pr(>Chisq)`[2]

  lmm_adj = lme4::lmer(Nadj ~ age + coverage + (1 | donor),
                       data = df_mut_count)
  lmm_null_adj = lme4::lmer(Nadj ~ age + (1 | donor), data = df_mut_count)
  anv_adj = anova(lmm_null_adj, lmm_adj) # Likelihood ratio test
  df_results[, 2] = anv_adj$`Pr(>Chisq)`[2]

  library(ggplot2)
  unadj_title = paste("1. unadjusted muts\neffect coverage p = ",
                      format(anv$`Pr(>Chisq)`[2], digits = 3))
  adj_title = paste("2. adjusted muts\neffect coverage p = ",
                      format(anv_adj$`Pr(>Chisq)`[2], digits = 3))

  df_mut_long = df_mut_count |>
      pivot_longer(cols = c(N, Nadj), values_to = "mutations", names_to = "adjustment") |>
      mutate(adjustment = ifelse(adjustment == "N", unadj_title, adj_title),
             donor = as.factor(donor))
  plot = df_mut_long |>
      ggplot(aes(x = coverage, y = mutations, color = donor)) +
      geom_point(alpha = 0.5) +
      stat_smooth(geom = "line", method = "lm", alpha = 0.5) +
      facet_grid(. ~ adjustment) +
      theme_bw() +
      scale_y_continuous(limits = c(0, NA), labels = scales::comma) +
      theme(legend.position = "none")

  return(list(df_results = df_results, plot = plot))
}
