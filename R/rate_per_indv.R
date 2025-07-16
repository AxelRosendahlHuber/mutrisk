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


# return the mutation rate for each individual sample, splitted out by signature
sig_rate_per_indv = function(single_mut_sig_rates, sig_contribution) {

  mutrates = single_mut_sig_rates |>
    left_join(triplet_match_substmodel, by = c("triplet", "strand")) |>
    select(mut_type, signature, mle, cilow, cihigh)

  sig_contribution = sig_contribution |>
    column_to_rownames("sampleID") |>
    as.matrix()

  cell_rate_list = list()
  for (col in c("mle", "cilow", "cihigh")) {
    print(col)

    mutrates_col = mutrates |>
      select(signature, mut_type, all_of(col)) |>
      pivot_wider(values_from = all_of(col), names_from = signature) |>
      column_to_rownames("mut_type")

    # calculate the trinucleotide rate for each of the patients in the cohort
    sig_rates = lapply(as.data.frame(t(sig_contribution)), \(x) as.data.table(x* t(mutrates_col), keep.rownames = "signature"))
    cell_rate_list[[col]] = rbindlist(sig_rates, idcol = "sampleID") |>
      pivot_longer(-c(signature, sampleID), names_to = "mut_type", values_to = col)
  }

  # merge the low and high confidence intervals
  cell_rate_list[[1]]$cilow = cell_rate_list[[2]]$cilow
  cell_rate_list[[1]]$cihigh = cell_rate_list[[3]]$cihigh

  return(cell_rate_list[[1]])
}

