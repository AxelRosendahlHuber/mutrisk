#' Calculate sample specific mutation rates
#'
#' @param single_mut_sig_rates Mutation rates, output of \code{\link{get_mutrisk}}
#' @param sig_contribution Contribution of different signatures to the mutations of samples.
#' Output of \code{\link{fit_muts_to_sigs}}
#'
#' @returns A data.frame containing the mutation probabilities for a single cell in a given individual.
#'  with the columns:
#'  - sampleID: the name of the sample
#'  - mut_type: mutation type for which the specific probability is taken
#'  - mle: maximum likelihood estimate, the value which is typically used
#'  - cilow: 95% lower confidence interval according to the mutation rate model
#'  - cihigh: 95% higher confidence interval according to the mutation rate model
#' @export
#'
#' @examples
rate_per_indv = function(single_mut_sig_rates, sig_contribution) {

  mutrates = single_mut_sig_rates |>
    left_join(triplet_match_substmodel, by = c("triplet", "strand")) |>
    select(mut_type, signature, mle, cilow, cihigh)

  sig_contribution = sig_contribution |>
    column_to_rownames("sampleID") |>
    as.matrix()

  cell_rate_list = list()
  for (col in c("mle", "cilow", "cihigh")) {

    mutrates_col = mutrates |>
      select(signature, mut_type, all_of(col)) |>
      pivot_wider(values_from = all_of(col), names_from = signature) |>
      column_to_rownames("mut_type")

    # calculate the trinucleotide rate for each of the patients in the cohort
    cell_rate_list[[col]] = apply(sig_contribution, 1, \(x) colSums(t(mutrates_col) * x))  |>
      as.data.table(keep.rownames = "mut_type") |>
      pivot_longer(-mut_type, names_to = "sampleID", values_to = col)
  }

  # return the shared mutation rate list
  left_join(cell_rate_list[[1]], cell_rate_list[[2]], by = c("sampleID", "mut_type")) |>
    left_join(cell_rate_list[[3]], by = c("sampleID", "mut_type")) |>
    dplyr::select(sampleID, mut_type, mle, cilow, cihigh)
}


#' Calculate sample specific mutation rates, specific by signature
#'
#' @param single_mut_sig_rates Mutation rates, output of \code{\link{get_mutrisk}}
#' @param sig_contribution Contribution of different signatures to the mutations of samples.
#' Output of \code{\link{fit_muts_to_sigs}}
#'
#' @returns A data.frame containing the mutation probabilities for a single cell in a given sample.
#'  with the columns:
#'  - sampleID: the name of the sample
#'  - mut_type: mutation type for which the specific probability is taken
#'  - signature: the active mutational signature
#'  - mle: maximum likelihood estimate, the value which is typically used
#'  - cilow: 95% lower confidence interval according to the mutation rate model
#'  - cihigh: 95% higher confidence interval according to the mutation rate model
#' @export
#'
#' @examples
sig_rate_per_indv = function(single_mut_sig_rates, sig_contribution) {

  mutrates = single_mut_sig_rates |>
    left_join(triplet_match_substmodel, by = c("triplet", "strand")) |>
    select(mut_type, signature, mle, cilow, cihigh)

  sig_contribution = sig_contribution |>
    column_to_rownames("sampleID") |>
    as.matrix()

  mutrates_mle = mutrates |>
      select(signature, mut_type,mle) |>
      pivot_wider(values_from = mle, names_from = signature) |>
      column_to_rownames("mut_type")

  # calculate the trinucleotide rate for each of the patients in the cohort
  sig_rates = lapply(as.data.frame(t(sig_contribution)), \(x) as.data.table(x* t(mutrates_mle), keep.rownames = "signature"))
  cell_rates = rbindlist(sig_rates, idcol = "sampleID")
  return(cell_rates)
}
