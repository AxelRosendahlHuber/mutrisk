#  
#' Cosine similarity
#' @description function copied from the R Bioconductor package "MutationalPatterns"
#' \link[MutationalPatterns:cos_sim]{cos_sim} function
#'  
#' @param x 
#' @param y 
#'
#' @returns
#' @export
#'
#' @examples
cos_sim = function(x, y) {
  res <- x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  res <- as.numeric(res)
  return(res)
}
