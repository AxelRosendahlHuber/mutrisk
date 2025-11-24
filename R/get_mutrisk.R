#' Get signature-specific mutation probabilities
#'
#' @param muts_context list of input mutations with context provided, output of \code{\link{get_mut_context}}
#' @param dnds dNdS output list (can be the output of either the `dndscv` (dndscv package) or
#'  `dndscv_intron` (wintr package) functions)
#' @param sig_contribution data.frame with the signature contribution with sampleID column indicating the sample,
#' @param input_signatures Matrix or dataframe with signatures
#'  and other columns being numeric columns indicating the active mutational sigantures
#' @returns
#' @export
#'
#' @examples
get_mutrisk = function(muts_context, dnds,  sig_contribution, input_signatures) {

  # convert signature contribution data.frame to matrix
  if ("sampleID" %in% colnames(sig_contribution)) {
    sig_contribution = sig_contribution |>
      column_to_rownames("sampleID")
  }


  # compare exome and whole-genome mutation data
  # quick check if the ratio observed/annotated mutations is not too high - warning that this could be a sign of using exome data
  ratio_obs_annotated = nrow(muts_context)/nrow(dnds$annotmuts |> filter(gene != "intronic")) < 10
  if (ratio_obs_annotated) {
    print(paste0("ratio observed mutations/annotated mutations is: ", round(ratio_obs_annotated,1),
                 "\n possible due input data being exome or targeted data"))
  }

  # get the rate across the cohort for each signature (all activities summed up together)
  rates_by_sig = muts_to_sig(dnds, sig_contribution, input_signatures)

  # compare rates of annotated mutations in RefCDS and all mutations. To 'simulate' mutations in a given example using the input data, a correction is necessary
  sig_all_rates = data.frame(wgs_sig_activity = colSums(sig_contribution),
                             signature = colnames(sig_contribution))

  # get the trinucleotide counts for the observed regions:
  # the mutation rates in rates_by_sig are the total mutation rates observed in the cohort in the refCDS, calculated by trinucleotide
  # 1: The mutation rate for the trinuc level in refcds = rates_by_sig
  # 2: The mutation rate for the whole-genome level
  signatures = input_signatures[, sig_all_rates$signature] # pre-select active signatures present in the tissue
  wg_trinuc_rates = t(t(signatures) *  sig_all_rates$wgs_sig_activity) |>
    as.data.frame() |>
    rownames_to_column("triplet")  |>
    left_join(triplet_match_substmodel |> select(trinuc, triplet) |> distinct(), by = "triplet") |>
    left_join(hg19_trinuc_counts |> dplyr::rename(trinuc = trinucleotide), by = "trinuc") |>
    mutate(across(starts_with("SBS"), ~ . / trinuc_counts))


  # 3. The mutation rate for the intergenic level
  # trinuc counts whole genome
  # short part calculating three levels of different counts:
  # rates for a single signature
  refcds_trinuc_counts = trinuc_counts_dnds(dnds$L) # count the tricnucleotide counts present in the exonic regions

  # correct the input rates for the all rates (wgs data: whole-genome, targeted data: all regions) rates and calculate the mutation rate for a single mutation
  # mutation rates are divided by the total signature contribution
  single_mut_sig_rates = left_join(rates_by_sig, sig_all_rates, by = "signature") |>
    mutate(across(c(mle, cilow, cihigh), \(x) x / wgs_sig_activity)) |>
    select(signature, triplet, strand, mle, cilow, cihigh)

  #### Calculate the mean across the individual cells across patients
  cell_contribution = sig_contribution |>
    rownames_to_column("sampleID") |>
    pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "wgs_muts_cell")

  single_cell_rates = single_mut_sig_rates |>
    left_join(refcds_trinuc_counts, by = c("triplet", "strand")) |>
    left_join(cell_contribution, relationship = "many-to-many", by = "signature") |>
    mutate(mle = mle * trinuc_counts * wgs_muts_cell) |>
    group_by(sampleID, signature, wgs_muts_cell) |>
    summarize(exonic_muts_cell = sum(mle), .groups = "drop")

  mutrisk_results = list(single_mut_sig_rates = single_mut_sig_rates,
                         single_cell_rates = single_cell_rates)

  return(mutrisk_results)
}
