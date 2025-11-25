#' min_entropy function
#' 
#' Computes the minimum entropy score of an alignment matrix. 
#'
#' @param A : Character matrix for an alignment
#' @param pseudocount : Numeric scalar. Pseudocount for estimating the probabality for every letter at every position. Default = 1.
#'
#' @return : Minimum entropy score for alignment.
#' @export
#'
#' @examples A = matrix(c("A", "A", "T", "C", "A","A", "C", "C", "T", 
#' "G","A", "A", "-", "C", "T"), nrow = 3, byrow = TRUE)
#' min_entropy(A)
#' 
#' 
source("./functions/get_column_stats.R")

min_entropy <- function(A, pseudocount = 1){
  
  #gets counts (c_ia's) and probs (p_ia's)
  res = get_column_stats(A, pseudocount = pseudocount)
  C = res$count_mat
  P = res$p_mat
  
  nsites = ncol(C)
  
  scores = rep(NA, nsites)
  #loop through columns and compute score
  for(i in 1:nsites){
    c_vec = C[,i]
    p_vec = P[,i]
    colscore = -c_vec*log(p_vec)
    scores[i] = sum(colscore)
  }
  
  total_score = sum(scores)
  return(total_score)
}

testthat::test_that("",{
  A = matrix(
    c("A", "A", "B",
      "A", "B", "C",
      "A", "A", "B"),
    nrow = 3, byrow = TRUE)
  
  output = min_entropy(A, pseudocount = 1)
  
  #get scores for every colunm manually
  ent = rep(NA,3)
  ent[1] = 3*log(4/6)
  ent[2] = 2*log(3/6) + 1*log(2/6)
  ent[3] = 2*log(3/6) + 1*log(2/6)
  ans = sum(ent)*-1
  
  testthat::expect_equal(output, ans)
})

testthat::test_that("",{
  A = matrix(
    c("A", "A", "B","-",
      "A", "B", "C","B",
      "A", "A", "B","B"),
    nrow = 3, byrow = TRUE)
  
  output = min_entropy(A, pseudocount = 1)
  
  #get scores for every colunm manually
  ent = rep(NA,4)
  ent[1] = 3*log(4/7)
  ent[2] = 2*log(3/7) + 1*log(2/7)
  ent[3] = 2*log(3/7) + 1*log(2/7)
  ent[4] = 1*log(2/7) + 2*log(3/7)
  ans = sum(ent)*-1
  
  testthat::expect_equal(output, ans)
})
