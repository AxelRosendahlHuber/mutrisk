#' #' Pipeline modeling mutation rates and performing side-analyses
#' #'
#' #' Pipeline function chaining the functions of the mutrisk package
#' #'  1. Estimate the neutral mutation rates across the cohort
#' #'  2. Performing
#' #'
#' #'
#' #' @param cell_muts Dataframe containing mutations, in the format: sampleID (ID of individual sample, can be more for each donor)
#' #' chr (chromosome), pos (genomic position), #' ref (reference allele), alt (mutated allele), category (in case there are multiple different categories), donor (patient/donor ID)
#' #' @param select_sigs Optional: select specific signatures from the COSMIC reference signaturesv3.4: https://cancer.sanger.ac.uk/signatures/sbs/ [default: none]
#' #' @param tissue Indicate the tissue being studied
#' #' @param WES Logical, set to TRUE if using whole-exome data [default: FALSE]
#' #' @param multiple_refit_methods Logical - is signature analysis required for multiple methods? [default: FALSE]. For methods being
#' #' tested in the extensive analysis see \code{\link{fit_all_sig_methods}}
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' process_mutations = function(cell_muts, select_sigs, tissue, WES = FALSE, multiple_refit_methods = FALSE) {
#'
#'   # [inputs]
#'   outdir = paste0("processed_data/", tissue, "/")
#'
#'   if (WES) {
#'     region = "WES"
#'   } else {region = "WGS"}
#'
#'   names = c("sampleID", "chr", "pos", "ref", "alt", "category", "donor")
#'   if (!all(names %in% colnames(cell_muts))) {
#'     stop(paste("colnames muts need to be", paste(names, collapse = ", ")))
#'   }
#'
#'   # if provided, select specific signatures
#'   if (!missing(select_sigs)) {
#'     signatures = signatures[,select_sigs] # pre-select active signatures present in the tissue
#'   }
#'
#'   muts = as.data.frame(cell_muts)[,names]
#'
#'   # 1. perform dNdScv and determine the mutation rates of the individual processes:
#'
#'   message("1. model mutation rates using dNdScv and dNdScv-intron (wintr)")
#'   list_muts = split(muts, muts$category)
#'   annotmuts_list = list()
#'   sig_contribution_list = list()
#'   mut_dist_list = list()
#'   mut_dist_type_list = list()
#'
#'   for (select_category in names(list_muts)) {
#'     print(select_category)
#'     category_muts = list_muts[[select_category]]
#'     unique_muts = category_muts |> dplyr::select(-sampleID) |>
#'       dplyr::rename(sampleID = donor) |> distinct() |>   # select unique mutations by donor to exclude shared mutations across cell lineages.
#'       select(sampleID, chr, pos, ref, alt)
#'     dnds_intron_unique = dndscv_intron(unique_muts, outmats = TRUE, refdb = region, transcript_regions = region)
#'     save_dnds(dnds_intron_unique, folder = paste0(outdir,"/unique_", tissue, "_", select_category, "_"), outmats = TRUE, selection = TRUE)
#'
#'     # perform extra dNdS analysis using all mutations (not only unique)
#'     # to annotate individual cells, shared mutations can be added.
#'     dnds_intron = dndscv_intron(list_muts[[select_category]], outmats = TRUE,  refdb = region, transcript_regions = region)
#'     save_dnds(dnds_intron, folder =  paste0(outdir,"/", tissue, "_", select_category, "_"), outmats = TRUE)
#'
#'     # perform the 'regular' dnds to get the mutation rates when using exonic muations only
#'     dnds_exon = dndscv(unique_muts,  outmats = TRUE)
#'     save_dnds(dnds_exon, folder =  paste0(outdir,"/exon_", tissue, "_", select_category, "_"), outmats = TRUE)
#'
#'     ### signature extraction #####
#'     # perform signature extraction analysis on all mutations
#'     message("2. Refit mutational signatures to the existing mutational profiles")
#'     muts_context = get_mut_context(list_muts[[select_category]] |> dplyr::select(-donor))
#'     mm = muts_to_mat(muts_context)
#'
#'     if (multiple_refit_methods == TRUE) {
#'       message("setting for extensive signature analysis on: refitting using MuSiCal, MutationalPatterns and SigLasso")
#'       signature_refits = fit_all_sig_methods(mm, signatures)
#'       # plot and save the reconstruction of the mutational signatures for the individual samples
#'       plot_refits = plot_sig_extraction_methods(signature_refits)
#'       cossim_refits = cos_sim_reconstructed_methods(signature_refits, signatures, mm)
#'       contribution_cosine = compare_contribution_cosine(signature_refits)
#'       refit_check_plots = plot_refits / (cossim_refits | contribution_cosine | plot_spacer() ) +
#'         plot_layout(heights  = c(1.5, 1)) +
#'         plot_annotation(tag_levels = "A") &
#'         theme(plot.margin = margin(5,5,5,5, unit = "mm"),
#'               plot.tag = element_text(face = 'bold'))
#'
#'       ggsave(paste0("plots/", tissue, "/", select_category, "_", "refit_check.png"), refit_check_plots, width = 18, height = 14, bg = "white")
#'
#'       # add the signature extraction from the selected category as output. Take the MuSical refitted values as final value
#'       sig_contribution_list[[select_category]] = signature_refits |>
#'         filter(method == "MuSiCal") |>
#'         dplyr::select(-c(method, unassigned, total))
#'     } else {
#'         message("setting for extensive signature analysis off: refitting signatures using MuSiCal")
#'        sig_contribution_list[[select_category]] = fit_muts_to_sigs(mm, signatures, method = "MuSiCal")
#'     }
#'
#'     # Check the distribution of mutations across the genome
#'     mut_dist = check_mutation_distribution(mm, dnds_intron)
#'     ggsave(paste0("plots/", tissue, "/", select_category, "_", "mutation_distribution.png"), mut_dist$plot, width = 18, height = 14, bg = "white")
#'
#'     mut_dist_list[[select_category]] = mut_dist$mutrates_overall
#'     mut_dist_type_list[[select_category]] = mut_dist$mutrates_type
#'   }
#'
#'   # save signature contribution
#'   sig_contribution = rbindlist(sig_contribution_list)
#'   fwrite(sig_contribution, paste0(outdir, "/signature_contributions.tsv"), sep = "\t")
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
#'     })  |>
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
