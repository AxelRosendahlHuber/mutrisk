# function to get the trinucleotide counts from DNA sequences:
#' get trinucleotide sequences
#'
#' @param sequences DNAstringset to count the trinucleotide sequences(output of getSeq from a genome object)
#'
#' @returns dataframe with trinucleotides and their respective counts
#' @export
#'
#' @examples
get_trinuc_counts = function(sequences) {
  chrom_counts = trinucleotideFrequency(sequences)

  sums = colSums(chrom_counts)
  purine_index = substr(names(sums), 2,2) %in% c("A", "G")
  purines = sums[purine_index]
  pyrimidines = sums[!purine_index]

  names(purines) = names(purines) |>
    DNAStringSet() |>
    reverseComplement() |>
    as.character()

  trinuc_counts = purines[names(pyrimidines)] + pyrimidines
  as.data.frame(trinuc_counts) |>
    tibble::rownames_to_column("trinuc")
}
