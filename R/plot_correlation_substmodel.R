# plot correlation rate for 1, 12, 96 and 192 trinucleotide rates across the individual sampesl

plot_correlation_substmodel = function(dnds_intron) {

  # get the average mutation rate compared against trinucs:
  df_syn_rates = triplet_match_substmodel
  flat_rate = sum(dnds_intron$N)/sum(dnds_intron$L)

  # observed_rates- (mean of intron + syn)

  N = rowSums(dnds_intron$N, dims = 2)
  L = rowSums(dnds_intron$L, dims = 2)
  muts_syn = N[,1]
  muts_intr = N[,5] / dnds_intron$globaldnds[6,2]
  obs_rates = (muts_syn + muts_intr) / L[,5]

  # 192 submodel rates
  rates_192 = get_syn_mutrates(dnds_intron)

  # 12 base mutation rates
  submodel = dnds_intron$mle_submodel_12
  rates =  setNames(submodel[,2], submodel[,1])
  rates_12 = sapply(substmodel_introns_12[,1], function(x) prod(rates[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site

  rate_model = df_syn_rates |>
    mutate(obs_rates = obs_rates,
             flat_rate = flat_rate,
           rates_192 = rates_192,
           rates_12 = rates_12) |>
    group_by(triplet) |>  mutate(rate_96 = mean(rates_192)) |>
    group_by(type) |>  mutate(rate_6 = mean(rates_12))  |>
    ungroup() |>
    select(-strand, -trinuc, -mut, triplet, -type, -triplet)  |>
    column_to_rownames("mut_type")


  rate_model_diff = rate_model |>
    mutate(across(everything(),  ~ . / obs_rates)) |>
    rownames_to_column("mut_type") |>
    left_join(triplet_match_substmodel) |>
    dplyr::select(-obs_rates)

  rate_model_diff |>
    pivot_longer(c(rates_192, rate_96, rates_12, rate_6, flat_rate)) |>
    mutate(name = factor(name, levels = c("flat_rate", "rate_6", "rates_12", "rate_96", "rates_192"))) |>
    ggplot(aes(x = name, y = value, color = name)) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom() +
    geom_hline(yintercept = 1, linetype = "dotted") +
    facet_wrap(. ~ type) +
    theme_classic() +
    cowplot::panel_border() +
    ggsci::scale_color_igv() +
    labs(y = "error towards average observed 192-trinucleotide rate", x = NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}



# todo: update the function with two additional capabilities:
# 1: correct for the relative mutagenicity (intergenic, whole genome)
# 2: correct for the relative mutagenicity rate for a specific gene
# 3. Check if possible to do this for a list of mutations, in order to model the mutation rate for a specific gene.


