#' Convert list of mutation contexts to 96-trinucleotide matrix
#'
#' @param muts_context Data frame, data table or tibble with mutations ordered by sampleID, triplet (96 categories)
#'
#' @returns ordered matrix containing the rows the 96-trinucleotides and the columns the different samples present in the sampleID column
#' @export
#'
#' @examples
muts_to_mat = function(muts_context) {

  mut_mat = muts_context |>
    group_by(sampleID, triplet) |>
    count() |>
    pivot_wider(names_from = "triplet", values_from = n, values_fill = 0)
  mm = mut_mat |>
    column_to_rownames("sampleID")  |>
    as.matrix() |> t()

  out_mm = matrix(0, ncol = ncol(mm), nrow = 96, dimnames = list(TRIPLETS_96, colnames(mm)))
  out_mm[rownames(mm),] = mm

  return(out_mm)
}
