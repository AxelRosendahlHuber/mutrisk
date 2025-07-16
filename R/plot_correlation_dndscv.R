plot_correlation_dnsdcv = function(dnds_intron, L) {
  
  avg_rate = sum(dnds_intron$N)/sum(dnds_intron$L)
  gm = dnds_intron$genemuts # shorthand for genemuts
  ratio_cv = gm$exp_syn_cv / gm$exp_syn
  
  df_syn = data.table(
    observed_synonymous = gm$n_syn,
    length = colSums(L[,1,]) * avg_rate,
    `length + covs` = colSums(L[,1,]) * avg_rate * ratio_cv,
    `length + trinuc` = gm$exp_syn,
    `length + trinuc + covs` = gm$exp_syn_cv)
  
    data.frame(correlation = cor(df_syn)[-1,1]) |>  
      mutate(name = colnames(df_syn)[-1]) |> 
      ggplot(aes(x = name, y = correlation)) + 
      geom_col() + 
      theme_classic() + 
      scale_y_continuous(expand=expansion(mult=c(0,0.1))) + 
      geom_text(aes(label = format(round(correlation, 2), trim = T, drop0trailing = TRUE, big.mark = ","), 
                    y = correlation), vjust = -0.2, position = position_dodge(0.9)) + 
      labs(x = NULL)
}


plot_correlation_dnsdcv_intron = function(dnds_intron, L) {
  
  avg_rate = sum(dnds_intron$N)/sum(dnds_intron$L)
  gm = dnds_intron$genemuts # shorthand for genemuts
  ratio_cv = gm$exp_syn_cv / gm$exp_syn
  
  wintr = dnds_intron$globaldnds[6,2]
  
  length = (colSums(L[,1,]) + (colSums(L[,5,]) *wintr))
  df_syn = data.table(
    observed_synonymous = gm$n_syn + gm$n_intr,
    length = length * avg_rate,
    `length + covs` = length * avg_rate * ratio_cv,
    `length + trinuc` = gm$exp_syn + gm$exp_intr,
    `length + trinuc + covs` = gm$exp_intr_syn_cv)
  
  data.frame(correlation = cor(df_syn)[-1,1]) |>  
    mutate(name = colnames(df_syn)[-1]) |> 
    ggplot(aes(x = name, y = correlation)) + 
    geom_col() + 
    theme_classic() + 
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) + 
    geom_text(aes(label = format(round(correlation, 2), trim = T, drop0trailing = TRUE, big.mark = ","), 
                  y = correlation), vjust = -0.2, position = position_dodge(0.9)) + 
    labs(x = NULL)
}
