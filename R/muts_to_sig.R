#' Split mutation rates by signature
#'
#' @param dnds dnds object returned from \link[wintr::dndscv_intron]{dndscv_intron}.
#' For exome, output from \link[dndscv::dndscv]{dndscv} can be used as well
#' @param sig_contribution data frame of signature contributions. rownames indicating samples, columns indicating the signature activity
#' (either relative or absolute values)
#' followed by columns indicating the values for different signatures, with each column denoting a single signature
#' @param input_signatures Optional: Provide different input signatures [default: Cosmic v3.4 signatures]
#'
#' @returns Mutation rates split by signature activity
#' @export
#'
#' @examples
muts_to_sig = function(dnds, sig_contribution, input_signatures) {

  if (!missing(input_signatures)) {
    signatures = input_signatures # use the user-provided signatures
    signatures = signatures[TRIPLETS_96,] # make sure the signatures are formatted in the right way
  }

  if (!(all(colnames(sig_contribution) %in% colnames(signatures)))) {
    stop("signature names of signature contribution and profile matrix do not match")
  }

  rel_total_sigs = colSums(sig_contribution) |> proportions() # normalize the signature contributions
  sigs_select = signatures[, names(rel_total_sigs)] # select only the signatures used for the cohort

  # multiply the profiles of the signatures to their respective activities in the cohort.
  # Then, normalize along trinucleotides to get the relative contribution of each signature to a trinucleotide
  sigs_prop = t(t(sigs_select) * rel_total_sigs)
  sigs_rowfract = prop.table(sigs_prop, 1) |>
    as.data.table(keep.rownames = "triplet")

  prediction_list = vector("list", length = 3)
  names(prediction_list) = c("mle", "cilow", "cihigh")

  index_zero = rowSums(dnds$N) == 0
  trinucs_zero = mut_types[index_zero]
  if (sum(index_zero) > 0) {
    message(paste0(sum(index_zero),
    "trinucleotides with 0 mutations - impossible to calculate trinucleotide-specific mutation rates
    This rate of the following trinucleotides will be replaced by the rate of a 12-type substitution model\n",
    paste(trinucs_zero, collapse = "\n")))
  }

  # define function th calculate the mutation rates based on the substitution model
  calculate_mutrates <- function(mle_submodel, substmodel) {
    parmle <- setNames(mle_submodel[, col], mle_submodel[, 1])
    sapply(substmodel[, 1], \(x) prod(parmle[strsplit(x, "\\*")[[1]]]))
  }

  for (col in names(prediction_list)) {

    mutrates = calculate_mutrates(dnds$mle_submodel, substmodel)
    mutrates_12 = calculate_mutrates(dnds$mle_submodel_12, wintr::substmodel_introns_12)
    mutrates[index_zero] = mutrates_12[index_zero]  # replace all trinucleotides for which no mutations has been found to the 12-substitution model

    mutrates_df = data.frame(mut_type = names(mutrates), mutrate = mutrates)
    mutrate_match = left_join(mutrates_df, triplet_match_substmodel, by = "mut_type")

    # merge the rates, and multiply each signature column with the mutrate to get signature-specific mutation rates
    sig_mutrate = left_join(mutrate_match, sigs_rowfract, by = "triplet") |>
      mutate(across(starts_with("SBS"), \(x) x * mutrate)) |>
      dplyr::select(-trinuc)

    prediction_list[[col]] = sig_mutrate
    signatures = colnames(sig_mutrate |> dplyr::select(starts_with("SBS")))

  }

  rates_by_sig = rbindlist(prediction_list, idcol = "pred_type") |>
    select(-mutrate) |>
    pivot_longer(starts_with("SBS"), names_to = "signature") |>
    pivot_wider(names_from = "pred_type", values_from = "value") # get the signature scores

  return(rates_by_sig)
}
