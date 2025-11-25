#' consensus_score function
#' 
#' Takes a proportions matrix (output of get_column_proportions) and computes consensus score. 
#' Score = average consensus over all the sites/number of sequences
#'
#' @param prop_mat : List(Numeric matrix of proportions for letters for every column in an alignment matrix, number of rows in alignment matrix). Output of get_column_proportions().
#'
#' @return : Numeric vector of raw consensus values for every column in prop_mat. 
#' @export
#'
#' @examples see example

source("./functions/get_column_proportions.R")
consensus_score <- function(prop_mat){

  #get consensus of every column
  M = prop_mat$count_mat
  cols = split(M, rep(1:ncol(M), each = nrow(M)))
  con = sapply(cols, max)
  return(con)
}

#test ----

A = matrix(
  c("A", "A", "T", "C", "A", "T", "A",
    "C", "C", "C", "T", "G", "-", "A",
    "A", "A", "-", "C", "T", "-", "A"),
  nrow = 3, byrow = TRUE
)

X = get_column_proportions(A)
output = consensus_score(X)

C = c(2/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1)
#test the raw scores
output = consensus_score(X)
testthat::expect_equal(as.numeric(output), C)

