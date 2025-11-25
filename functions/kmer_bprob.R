#' kmer_bprob
#' 
#' Compute the probability of generating a kmer under the random model (with background probabilities)
#'
#' @param kmer : A string vector of kmer
#' @param background : NumericMatrix of background letter probabilites
#' @param log_prob: Logical. Set to TRUE to return log probabilities. Default to F. 
#'
#' @return Numeric scalar. Probability of kmer under background model. 
#' @export
#'
#' @examples k = list(c("C","A", "B"))
#' q = c(0.2, 0.5, 0.3)
#' rownames(q) = c("A","B", "C")
kmer_bprob <- function(kmer, background, log_prob = F){
  #initialise vector of probs 
  prob = rep(NA, length(kmer))
  #loop through letters in the kmer and get probabilities
  for(pos in 1:length(kmer)){
    sym = kmer[pos]
    # If symbol not in profile columns, probability 0
    if(!(sym %in% rownames(background))){
      prob[pos] = 0
    } else {
      sym_index = which(rownames(background)==sym)
      prob[pos] = background[sym_index]
    }
  }
  #return log probs if log is T
  if(log_prob){
    lp = log(prob)
    kmer_prob = sum(lp)
  } else {
    #take the product of all positions
    kmer_prob = prod(prob)
  }
  return(kmer_prob)
}

#vanilla case
testthat::test_that("",{
  k = c("C","A", "B","B")
  q = matrix(c(0.1, 0.2, 0.7))
  rownames(q) = c("A","B","C")
  output = kmer_bprob(k, q)
  ans = 0.7*0.1*0.2^2
  testthat::expect_equal(output,ans)
})

#log case
testthat::test_that("",{
  k = c("C","A", "B","B")
  q = matrix(c(0.1, 0.2, 0.7))
  rownames(q) = c("A","B","C")
  output = kmer_bprob(k, q,log_prob = T)
  ans = log(0.7)+log(0.1)+2*log(0.2)
  testthat::expect_equal(output,ans)
})
