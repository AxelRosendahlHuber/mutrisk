# calculate the trinuc counts in the RefCDS (to correct for the actual rate)
trinuc_counts_dnds = function(dnds_array) {

  # "stack" all the trinucleotide frames
  all = abind::abind(dnds_array, along = 3)
  total = apply(all, c(1,2), sum)
  total = total[,1:4] # only select genic areas, exclude possible intronic parts

  # count the recurring trinucleotides
  trinuc_mat = as.data.frame(total) |>
    mutate(mut_type = rownames(substmodel),
           trinucleotide = substr(mut_type, 1,3)) |>
    arrange(trinucleotide)

  trinuc_mat$trinuc_counts = rowSums(total)
  trinuc_rates = left_join(trinuc_mat |> select(-starts_with("V")),
                           triplet_match_substmodel ,by = "mut_type")

  return(trinuc_rates)
}
