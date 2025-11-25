#' profiles function
#' 
#' Compute the probabilities of all letters in an alignment matrix.
#'
#' @param A : Character matrix for alignment. 
#' @param alphabet : Character vector indicating all possible letters.
#' @param pseudocount : Pseudocount for the computing probabilities. Default = 1. 
#'
#' @return A NumericMatrix showing the probabilities of all letters at every column of A. 
#' @export
#'
#' @examples A = rbind(c("A","B","B","A"), c("A","A","B","A"),c("A","C","B","C"))
#' alphabet = c("A","B","C")
#' profile(A, alphabet)

source("./functions/count_letters.R")

profile <- function(A, alphabet, pseudocount = 1){
  
  nsites = ncol(A)
  nletters = length(alphabet)
  nseq = nrow(A)
  
  #initialize profile matrix P
  Q = matrix(nrow = nletters, ncol = ncol(A))
  rownames(Q) = alphabet
  
  #loop through columns of A
  for(j in 1:nsites){
    #get counts of every letter in the column (site)
    cij = count_letters(letters = A[,j], alphabet = alphabet)
    #add entry to Q
    qis = (cij + pseudocount)/(nseq -1 + pseudocount*nletters)
    Q[,j] = qis
  }
  #return profile matrix Q
  return(Q)
}

testthat::test_that("",{
  #params
  A = rbind(c("A","B","B","A"), c("A","A","B","A"),c("A","C","B","C"))
  alphabet = c("A", "B", "C")
  pseudocount = 1
 
  #compute manual solution
  Q = matrix(0, ncol = ncol(A), nrow = length(alphabet))
  rownames(Q) = alphabet
  N = nrow(A)
  B = pseudocount*length(alphabet)
  
  #manual check
  Q[,1] = unlist( lapply(c(3,0,0), function(x){(x+pseudocount)/(N - 1 + B)}) )
  Q[,2] = unlist( lapply(c(1,1,1), function(x){(x+pseudocount)/(N - 1 + B)}) )
  Q[,3] = unlist( lapply(c(0,3,0), function(x){(x+pseudocount)/(N - 1 + B)}) )
  Q[,4] = unlist( lapply(c(2,0,1), function(x){(x+pseudocount)/(N - 1 + B)}) )
  
  output = profile(A, alphabet, pseudocount)
  
  testthat::expect_equal(Q, output)
})
