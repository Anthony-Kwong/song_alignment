source("./functions/count_letters.R")
#' get_column_stats function
#' 
#' Takes an character matrix and returns the counts and probabilities for all the letters for every column.
#' The counts are raw counts. The estimated probabilities can be adjusted with pseudocounts. Used for calculating
#' minimum entropy of an alignment. 
#'
#' @param sequence_matrix : A character matrix. 
#' @param pseudocount : Numeric scalar. Pseudocount. Default = 1.
#'
#' @return list(count_mat: matrix of counts, p_mat : matrix of proportions)
#' @export
#'
#' @examples See example tests. 
get_column_stats <- function(sequence_matrix, pseudocount = 1) {
  # Identify all unique letters across the matrix
  unique_letters <- sort(unique(as.vector(sequence_matrix)))
  nsites = ncol(sequence_matrix)
  nseqs = nrow(sequence_matrix)
  
  # Initialize an empty matrix to store counts
  counts_matrix <- matrix(0, nrow = length(unique_letters), ncol = nsites)
  rownames(counts_matrix) <- unique_letters
  #colnames(counts_matrix) <- paste0("Column", seq_len(ncol(sequence_matrix)))
  
  # Loop through columns and count every letter
  for (i in 1:nsites) {
    col_i <- sequence_matrix[,i]
    #count every letter among unique_letters
    raw_counts = count_letters(letters = col_i, alphabet = unique_letters)
    counts_matrix[,i] = raw_counts
  }
  
  #convert counts to proportions, accounting for pseudocount
  p_mat = (counts_matrix + pseudocount)/(nseqs + pseudocount*length(unique_letters))
  
  #tie both together as a list
  output = list(count_mat = counts_matrix, p_mat = p_mat)
  
  return(output)
}

testthat::test_that("",{
  # Test the function
  A = matrix(
    c("A", "A", "B", "C", "A",
      "A", "C", "C", "B", "D",
      "A", "A", "-", "C", "B"),
    nrow = 3, byrow = TRUE
  )
  
  output = get_column_stats(A, pseudocount = 1)
  
  #manual check on counts
  ans1 = matrix(0, nrow = 5, ncol = 5)
  ans1[,1] = c(0,3,0, 0, 0)
  ans1[,2] = c(0,2,0,1,0)
  ans1[,3] = c(1,0,1,1,0)
  ans1[,4] = c(0,0,1,2,0)
  ans1[,5] = c(0,1,1,0,1)
  rownames(ans1) = c("-","A","B","C","D")
  
  testthat::expect_equal(output$count_mat, ans1)
  
  ans2 = (output$count_mat + 1)/(nrow(A) + 5)
  testthat::expect_equal(output$p_mat, ans2)
})

