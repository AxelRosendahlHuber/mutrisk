#' #### EXPERIMENTAL ATTEMPT TO SIMPLIFY THE DATA
#' Check how to best organize this function
#'
#'
#' # uniform function to processes the mutations:
#'
#' #' R script combining the estimation of local mutation rates,
#' #'
#' #' TODO: Find a better name for the script - for now: pipeline_cohort_mutations
#' #'
#' #' @param cell_muts # dataframe containing mutations, in the format: sampleID (ID of individual sample, can be more for each donor)
#' #' chr (chromosome), pos (genomic position), #' ref (reference allele), alt (mutated allele), category (in case there are multiple different categories), donor (patient/donor ID)
#' #' @param select_sigs Optional: select specific signatures from the COSMIC reference signaturesv3.4: https://cancer.sanger.ac.uk/signatures/sbs/ [default: none]
#' #' @param WES Logical, set to TRUE if using whole-exome data [default: FALSE]
#' #' @param multiple_refit_methods Logical - is signature analysis required for multiple methods? [default: FALSE]. For methods being
#' #' tested in the extensive analysis see \code{\link{fit_all_sig_methods}}
#' #' @param outdir Output directory, typically named after the tissue of the analysis
#' #' @param prefix Optional prefix to specifiy the analysis. Can be used for a subtype of cells, or specific exposure
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#'
#' pipeline_cohort = function(cell_muts, select_sigs,  outdir, prefix, WES = FALSE, multiple_refit_methods = FALSE) {
#'
#'   # make output directory and path
#'   outdir_prefix = paste0(outdir, "/", prefix, "_")
#'
#'   if (WES) {
#'     region = "WES"
#'   } else {region = "WGS"}
#'
#'   names = c("sampleID", "chr", "pos", "ref", "alt", "donor")
#'   if (!all(names %in% colnames(cell_muts))) {
#'     stop(paste("colnames muts need to be", paste(names, collapse = ", ")))
#'   }
#'
#'   # if provided, select specific signatures
#'   if (!missing(select_sigs)) {
#'     signatures = signatures[,select_sigs] # pre-select active signatures present in the tissue
#'   }
#'
#'   output_list = list()
#'   plot_list = list()
#'
#'   muts = as.data.frame(cell_muts)[,names]
#'   # select unique mutations by donor to exclude shared mutations across cell lineages
#'   unique_muts = muts |> dplyr::select(-sampleID) |>
#'       dplyr::rename(sampleID = donor) |> distinct() |>
#'       select(sampleID, chr, pos, ref, alt)
#'   dnds_intron_unique = dndscv_intron(unique_muts, outmats = TRUE, refdb = region, transcript_regions = region)
#'   save_dnds(dnds_intron_unique, folder = paste0(outdir_prefix, "_unique_"), outmats = TRUE, selection = TRUE)
#'
#'   # perform extra dNdS analysis using all mutations (not only unique)
#'   # to annotate individual cells, shared mutations can be added.
#'   dnds_intron = dndscv_intron(muts, outmats = TRUE,  refdb = region, transcript_regions = region)
#'   save_dnds(dnds_intron, folder =  outdir_prefix, outmats = TRUE)
#'
#'   # perform the 'regular' dnds to get the mutation rates when using exonic muations only
#'   dnds_exon = dndscv(unique_muts,  outmats = TRUE)
#'   save_dnds(dnds_exon, folder =  paste0(outdir_prefix, "exon_"), outmats = TRUE)
#'
#'   ### signature extraction #####-----------------------------------------------------------------------
#'   # perform signature extraction analysis on all mutations
#'   muts_context = get_mut_context(muts |> dplyr::select(-donor))
#'   mm = muts_to_mat(muts_context)
#'
#'   if (multiple_refit_methods == TRUE) {
#'     message("setting for extensive signature analysis on: refitting using MuSiCal, MutationalPatterns and SigLasso")
#'     signature_refits = fit_all_sig_methods(mm, signatures)
#'     # plot the reconstruction of the mutational signatures for the individual samples
#'     plot_refits = plot_sig_extraction_methods(signature_refits)
#'     cossim_refits = cos_sim_reconstructed_methods(signature_refits, signatures, mm)
#'     contribution_cosine = compare_contribution_cosine(signature_refits)
#'     plot_list[["refit_check_plots"]] = plot_refits / (cossim_refits | contribution_cosine | plot_spacer() ) +
#'       plot_layout(heights  = c(1.5, 1)) +
#'       plot_annotation(tag_levels = "A") &
#'       theme(plot.margin = margin(5,5,5,5, unit = "mm"),
#'             plot.tag = element_text(face = 'bold'))
#'
#'   #ggsave(paste0(outdir_prefix, "_refit_check.png"), refit_check_plots, width = 18, height = 14, bg = "white")
#'
#'   # add the signature extraction from the selected category as output. Take the MuSical refitted values as final value
#'     output_list[["signature_refits"]] = signature_refits |>
#'         filter(method == "MuSiCal") |>
#'         dplyr::select(-c(method, unassigned, total))
#'
#'     } else {
#'       message("setting for extensive signature analysis off: refitting signatures using MuSiCal")
#'       output_list[["signature_refits"]] = fit_muts_to_sigs(mm, signatures, method = "MuSiCal")
#'     }
#'
#'     # Check the distribution of mutations across the genome
#'     mut_dist = check_mutation_distribution(mm, dnds_intron)
#'     plot_list[["mut_dist"]] = mut_dist$plot
#'     #ggsave(paste0("plots/", tissue, "/", select_category, "_", "mutation_distribution.png"), mut_dist$plot, width = 18, height = 14, bg = "white")
#'
#'     output_list[["mutrates_overall"]] = mut_dist$mutrates_overall
#'     output_list[["mutrates_type"]] = mut_dist$mutrates_type
#'   }
#'
#'   # # save signature contribution
#'   # sig_contribution = rbindlist(sig_contribution_list)
#'   # fwrite(sig_contribution, paste0(outdir, "/signature_contributions.tsv"), sep = "\t")
#'
#'   # save overview plots of the mutation distribution
#'   lapply(mut_dist_list, \(x) {x |> mutate(mutrate  = mutrate/ x$mutrate[3])})  |>
#'     rbindlist(idcol = "category") |>
#'     group_by(category) |>
#'     filter(region != 'whole_genome') |>
#'     ggplot(aes(x = region, y = mutrate, fill = category)) +
#'     geom_col(position = "dodge") +
#'     geom_text(aes(label = format(round(mutrate, digits = 2), big.mark = ","), y = mutrate),
#'               vjust = -0.2, position = position_dodge(0.9)) +
#'     labs(x = "region (whole genome = 1)", y = "relative mutrate to whole genome") +
#'     ggsci::scale_fill_aaas() +
#'     cowplot::theme_cowplot() +
#'     scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#'
#'   ggsave(paste0("plots/", tissue, "/mutrate_dist_all.png"), width = 6, height = 4.5, bg = "white")
#'
#'   # save overview plots of the mutation distribution
#'   mutation_dist_type = lapply(mut_dist_type_list, \(x) {x |>
#'       pivot_wider(values_from = mutrate, names_from = region) |>
#'       mutate(across(-type, ~ ./ whole_genome))
#'   })  |>
#'     rbindlist(idcol = "category")
#'
#'   fwrite(mutation_dist_type, paste0(outdir, "/mutrate_dist_region.tsv"), sep = "\t")
#'
#'   mutation_dist_type |>
#'     pivot_longer(cols = c("intergenic", "transcribed"),  names_to = "region", values_to = "mutrate") |>
#'     group_by(category) |>
#'     ggplot(aes(x = type, y = mutrate, fill = category, alpha = region)) +
#'     geom_col(position = "dodge") +
#'     facet_grid(category ~ . ) +
#'     cowplot::theme_cowplot() +
#'     scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#'     scale_alpha_manual(values = c(0.6, 1)) +
#'     labs(x = NULL, y = "relative_mutrate to whole genome")
#'   ggsave(paste0("plots/", tissue, "/mutrate_dist_type.png"), width = 10, height = 8, bg = "white")
#' }
#'
#' # combine multiple methods
#' # combine_plots
#'
#' # save plots function



