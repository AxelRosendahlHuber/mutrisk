# # #testing variables (remove once package is finished)
# input_signatures = c("SBS1", "SBS5", "SBS18", "SBS19")
# output_path = "~/Documents/"
# name = "mutrisk_test"
# WES = FALSE
# multiple_refit_methods = FALSE
# sensitivity_correction = TRUE
# cell_muts = mutrisk:::cell_muts
# triplet_match_substmodel = mutrisk:::triplet_match_substmodel
# library(tidyverse)
# library(wintr)
# source("data-raw/mutation_signature_variables.R")

#devtools::load_all()# load all the variables/functions


#' Pipeline modeling mutation rates and performing side-analyses
#'
#' Pipeline function chaining the functions of the mutrisk package
#'  1. Estimate the neutral mutation rates across the cohort using dNdScv and dNdScv-intron
#'  2. Signature extraction (optional to perform multiple refitting methods)
#'  3. Separate mutation rates by signature
#'  4. Combine them back by patient
#'
#' @param cell_muts Dataframe containing mutations, in the format: sampleID
#'  (ID of individual sample, can be more for each donor)
#' chr (chromosome), pos (genomic position), #' ref (reference allele), alt (mutated allele),
#'  category (in case there are multiple different categories), donor (patient/donor ID)
#' @param input_signatures Optional: select specific signatures from the
#'  COSMIC reference signatures v3.4: https://cancer.sanger.ac.uk/signatures/sbs/ [default: none]
#' @param output_path Path where output directory is made
#' @param metadata Optional: neccesary for obtaining patient-specific mutation rates.
#'  Contains at least 2 columns: sampleID and donor.
#' @param name Optional: Name of analysis (useful in case multiple categories/tissues are run).
#'  Will determine name of output folder and  files
#' @param WES Logical, set to TRUE if using whole-exome data [default: FALSE]
#' @param multiple_refit_methods Logical - is signature analysis required for multiple methods? [default: FALSE].
#'  For methods being tested in the extensive analysis see \code{\link{fit_all_sig_methods}}.
#' Note: refitting multiple methods for large datasets can be very compute-intensive
#' @param sample_sig_specific Is output per sample/signature desired (with many samples/signatures can return very large output)
#' @param sensitivity_correction Logical: is sensitivity analysis wanted.
#'  Requires also metadata with a 'sensitivity' column denoting the sensitivity.
#'
#' @returns
#' @export
#'
#' @examples
mutrisk_pipeline = function(cell_muts, input_signatures,
                            output_path, metadata = NULL,
                            name,  WES = FALSE,
                            multiple_refit_methods = FALSE,
                            sensitivity_correction = FALSE) {


  # Set output directory and make a subdirectory named "plots"
  if (exists("name")) {
    outdir = paste0(output_path, "/", name, "/")
  } else {
    outdir = paste0(output_path, "/mutrisk_analysis/")
  }

  if (!dir.exists(outdir)) { dir.create(outdir) }
  plotdir = paste0(outdir, "/plots/")
  if (!dir.exists(plotdir))  { dir.create(plotdir) }

  # Set output region if WES is set to TRUE
  if (WES) {
    region = "WES"
  } else {region = "WGS"}

  # Check names of input mutation file
  names = c("sampleID", "chr", "pos", "ref", "alt", "category", "donor")
  if (!all(names %in% colnames(cell_muts))) {
    stop(paste("colnames muts need to be", paste(names, collapse = ", ")))
  }
  muts = as.data.frame(cell_muts)[,names]

  # if provided, select specific signatures
  if (is.matrix(input_signatures)) {
    signatures = input_signatures
  }

  if (is.character(input_signatures)) {
    signatures = mutrisk::signatures # define the "signatures" variable as the active signatures in cosmic v3.4
    signatures = signatures[,input_signatures] # pre-select active signatures present in the tissue
  }

  # 1. perform dNdScv and determine the mutation rates of the individual processes:
  message("1. model mutation rates using dNdScv and dNdScv-intron (wintr)")
  category_muts = muts
  unique_muts = muts |> dplyr::select(-sampleID) |>
      dplyr::rename(sampleID = donor) |> dplyr::distinct() |>   # select unique mutations by donor to exclude shared mutations across cell lineages.
      dplyr::select(sampleID, chr, pos, ref, alt)
  dnds_intron_unique = wintr::dndscv_intron(unique_muts, outmats = TRUE, refdb = region, transcript_regions = region)
  wintr::save_dnds(dnds_intron_unique, folder = paste0(outdir,"/unique_", name, "_"), outmats = TRUE, selection = TRUE)

  # perform extra dNdS analysis using all mutations (not only unique)
  # to annotate individual cells, shared mutations can be added.
  dnds_intron = wintr::dndscv_intron(muts, outmats = TRUE,  refdb = region, transcript_regions = region)
  wintr::save_dnds(dnds_intron, folder =  paste0(outdir,"/", name, "_"), outmats = TRUE)

  # perform the 'regular' dnds to get the mutation rates when using exonic muations only
  dnds_exon = dndscv::dndscv(unique_muts,  outmats = TRUE)
  wintr::save_dnds(dnds_exon, folder =  paste0(outdir,"/exon_", name, "_"), outmats = TRUE)

  # inspect and plot correlations, intron/exon rates and correlations between observe and expected
  # mutation rates for genes
  if (region == "WGS") {
    L = abind::abind(lapply(wintr::RefCDS_WGS,\(x) x$L), along = 3)
  } else if (region == "WES") {
    L = abind::abind(lapply(wintr::RefCDS_WES,\(x) x$L), along = 3)
  } else { print("region argument must be either 'WES' or 'WGS'")}

  # plot correlations with mutation models
  plot_cor = plot_correlation_dnsdcv(dnds_intron, L) +
    labs(title = "Correlation expected synonymous to observed:", subtitle = name)
  ggsave(paste0(plotdir, name, "_corr_syn_obs.png"), plot_cor, width = 6, height = 3)

  plot_cor_intron = plot_correlation_dnsdcv_intron(dnds_intron, L) +
    labs(title = "Correlation expected synonymous + intron to observed:",subtitle = name)
  ggsave(paste0(plotdir, name, "_corr_syn_obs_intron.png"), plot_cor_intron, width = 6, height = 3)

  plot_cor_subs = plot_correlation_substmodel(dnds_intron) +
    labs(title = "Correlation expected synonymous to observed:", subtitle = name)
  ggsave(paste0(plotdir, name, "_substmodel.png"), plot_cor_subs, width = 10, height = 6)

  plot_comparison = compare_intron_exon(dnds_intron = dnds_intron, dnds_exon = dnds_exon) +
    plot_annotation(name) & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotdir, name, "_dnds_comparison.png"), plot_comparison, width = 16, height = 8)

  # get mutation matrix (also used in 2. signature Extraction)
  muts_context = get_mut_context(muts |> dplyr::select(-donor))
  mm = muts_to_mat(muts_context)

  # Check the distribution of mutations across the genome
  mut_dist = check_mutation_distribution(mm, dnds_intron)
  ggsave(paste0(plotdir, "/", name, "_mutation_distribution.png"), mut_dist$plot, width = 18, height = 14, bg = "white")
  fwrite(mut_dist$mutrates_overall, paste0(outdir, "mutrates_overall.txt"), sep = "\t")
  fwrite(mut_dist$mutrates_type, paste0(outdir, "mutrates_type.txt"), sep = "\t")

  ### signature extraction #####
  # perform signature extraction analysis on all mutations
  message("2. Refit mutational signatures to the existing mutational profiles")

  if (multiple_refit_methods == TRUE) {
    message("setting for extensive signature analysis on: refitting using
            1. MuSiCal,
            2. MutationalPatterns (standard)
            3. MutationalPatterns (strict)
            4. SigLasso")

    signature_refits = fit_all_sig_methods(mm, signatures)
    # plot and save the reconstruction of the mutational signatures for the individual samples
    plot_refits = plot_sig_extraction_methods(signature_refits)
    cossim_refits = cos_sim_reconstructed_methods(signature_refits, signatures, mm)
    contribution_cosine = compare_contribution_cosine(signature_refits)
    refit_check_plots = plot_refits / (cossim_refits | contribution_cosine | plot_spacer() ) +
      plot_layout(heights  = c(1.5, 1)) +
      plot_annotation(tag_levels = "A") &
      theme(plot.margin = margin(5,5,5,5, unit = "mm"),
            plot.tag = element_text(face = 'bold'))
      ggsave(paste0("plots/", tissue, "/", select_category, "_", "refit_check.png"), refit_check_plots, width = 18, height = 14, bg = "white")

      # add the signature extraction from the selected category as output. Take the MuSical refitted values as final value
      sig_contribution = signature_refits |>
        filter(method == "MuSiCal") |>
        dplyr::select(-c(method, unassigned, total))
    } else {
      message("setting for extensive signature analysis off: refitting signatures using MuSiCal")
      sig_contribution = fit_muts_to_sigs(mm, signatures, method = "MuSiCal")
    }

  # remove signatures for which no contribution is found
  sig_contribution = sig_contribution |>
      dplyr::select(sampleID, where(~ is.numeric(.) && any(. > 0)))
  # save signature contribution
  fwrite(sig_contribution, paste0(outdir, "/signature_contributions.tsv"), sep = "\t")

  message("3. Get the probabilities for each mutation of signature activity and for individual patients")
  mutrisk_rates = get_mutrisk(muts_context = muts_context, dnds = dnds_intron,
                              sig_contribution = sig_contribution)

  single_mut_sig_rates = mutrisk_rates$single_mut_sig_rates
  fwrite(single_mut_sig_rates, paste0(outdir, name, "_single_mut_sig_rates.tsv.gz"))

  # get the sample-specific single mutation rates
  rate_per_sample = rate_per_indv(single_mut_sig_rates, sig_contribution)
  # get the sample-specific signature mutation rates for each trinucleotide
  # (with many samples/signatures very large output)
  sig_rate_per_sample = sig_rate_per_indv(single_mut_sig_rates, sig_contribution)

  # 4. Optional: Correct for the relative sensitivity
  if (sensitivity_correction == TRUE & is.null(metadata)) {
    stop("missing metadata. Supply metadata with columns: sampleID, donor and sensitivity")
  }

  if (sensitivity_correction == TRUE & !is.null(metadata)) {
    message("4. Optional addition: Correct the mutation rates by sensitivity rate")

    if (!"sensitivity"  %in% colnames(metadata)) {
      stop("Could not correct the mutation rate values for sensitivity, as the  column 'sensitivity'
      is not present in the supplied metadata data.frame.
      Check the function get_sensitivity() to calculate the sensitivities across samples and run again")
    }

    sig_rate_per_sample = sig_rate_per_sample |>
      left_join(metadata, by = "sampleID") |>
      mutate(across(contains(">"), ~ ./sensitivity))

    rate_per_sample = rate_per_sample |>
      left_join(metadata, by = "sampleID") |>
      mutate(across(c(mle, cilow, cihigh), ~ ./sensitivity))

  }

  fwrite(rate_per_sample, paste0(outdir, name, "_rate_per_sample.tsv.gz"))
  fwrite(sig_rate_per_sample, paste0(outdir, name, "_sig_rate_per_sample.tsv.gz"))


  # output the data
  output_list = list(dnds_exon = dnds_exon,
                     dnds_intron_unique = dnds_intron_unique,  # consider not calculating this
                     dnds_intron = dnds_intron,
                     mutrisk_rates = mutrisk_rates,
                     rate_per_sample = rate_per_sample,
                     sig_rate_per_sample = sig_rate_per_sample)

  output_list
}

