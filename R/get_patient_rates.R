# # get the rates for individiual mutation rates.
# # also compute metrics for the influence of the mutation rate model
#
# pipeline_patient_rates = function(cell_muts, tissue, region = "WGS") {
#   outdir = paste0("processed_data/", tissue, "/")
#
#   categories = unique(cell_muts$category)
#   select_category = categories[1]
#
#   if (region == "WGS") {
#     L = abind::abind(lapply(wintr::RefCDS_WGS,\(x) x$L), along = 3)
#   } else if (region == "WES") {
#     L = abind::abind(lapply(wintr::RefCDS_WES,\(x) x$L), along = 3)
#   } else { print("region argument must be either 'WES' or 'WGS'")}
#
#   for (select_category in categories) {
#
#     print(select_category)
#     plotdir = paste0("plots/", tissue, "/", select_category, "/")
#     if (!(dir.exists(plotdir))) {dir.create(plotdir)}
#
#     dnds_intron_unique = readRDS(paste0(outdir ,"unique_", tissue, "_", select_category, "_dnds.rds"))
#     dnds_intron = readRDS(paste0(outdir, tissue, "_", select_category, "_dnds.rds"))
#     dnds_exon = readRDS(paste0(outdir, "exon_", tissue, "_", select_category, "_dnds.rds"))
#
#     # plot correlations with mutation models
#     plot_cor = plot_correlation_dnsdcv(dnds_intron, L) + labs(title = "Correlation expected synonymous to observed:",
#                                                            subtitle = paste0(tissue, "|", select_category))
#     ggsave(paste0("plots/", tissue, "/", select_category, "/", select_category, "_corr_syn_obs.png"), plot_cor, width = 6, height = 3)
#
#     plot_cor_intron = plot_correlation_dnsdcv_intron(dnds_intron, L) + labs(title = "Correlation expected synonymous + intron to observed:",
#                                                                          subtitle = paste0(tissue, "|", select_category))
#     ggsave(paste0("plots/", tissue, "/", select_category, "/", select_category, "_corr_syn_obs_intron.png"), plot_cor_intron, width = 6, height = 3)
#
#     plot_cor_subs = plot_correlation_substmodel(dnds_intron) + labs(title = "Correlation expected synonymous to observed:",
#                                                          subtitle = paste0(tissue, "|", select_category))
#     ggsave(paste0("plots/", tissue, "/", select_category, "/", select_category, "_substmodel.png"), plot_cor_subs, width = 10, height = 6)
#
#     plot_comparison = compare_intron_exon(dnds_intron = dnds_intron, dnds_exon = dnds_exon) +
#       plot_annotation(paste(tissue, "-", select_category)) & theme(plot.title = element_text(hjust = 0.5))
#     ggsave(paste0(plotdir, select_category, "_dnds_comparison.png"), plot_comparison, width = 16, height = 8)
#
#     muts = cell_muts |> filter(category == select_category)
#
#     sig_contribution = read.delim(paste0("processed_data/", tissue, "/signature_contributions.tsv")) |>
#       filter(sampleID %in% unique(muts$sampleID)) |>
#       column_to_rownames("sampleID")
#     sig_contribution = sig_contribution[,colSums(sig_contribution) > 0]
#
#     mutrisk_rates = get_mutrisk(muts = muts, dnds = dnds_intron,  metadata = metadata, sig_contribution = sig_contribution)
#
#     # get the individual single mutation rates
#     single_mut_sig_rates = mutrisk_rates$single_mut_sig_rates
#     fwrite(single_mut_sig_rates, paste0(outdir, select_category, "_single_mut_sig_rates.tsv.gz"))
#     rate_per_individual = rate_per_indv(single_mut_sig_rates, sig_contribution)
#     fwrite(rate_per_individual, paste0(outdir, select_category, "_patient_rates.tsv.gz"))
#
#     # get the signature-specific mutation rates for each trinucleotide
#     sig_rate_per_individual = sig_rate_per_indv(single_mut_sig_rates, sig_contribution)
#     fwrite(sig_rate_per_individual, paste0(outdir, select_category, "_sig_patient_rates.tsv.gz"))
#   }
# }
