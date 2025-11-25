# Helper: probability of a candidate k-mer under a profile
#' kmer_prob
#' 
#' Get probability of a candidate k-mer under a profile.
#'
#' @param kmer : A string vector of kmer
#' @param profile : NumericMatrix of profile of candidate alignment. 
#'
#' @return : Numeric. Probability of generating kmer under the given profile
#' @export
#'
#' @examples k = list(c("C","A", "B"))
#' q = cbind(c(0.25,0.5,0.25),c(0.5,0.5,0), c(0.25,0.5,0.25))
#' rownames(q) = c("A","B","C")

kmer_prob <- function(kmer, profile){
  #initialise vector of probs 
  prob = rep(NA, length(kmer))
  #loop through letters in the kmer and get probabilities
  for(pos in 1:length(kmer)){
    sym = kmer[pos]
    # If symbol not in profile columns, probability 0
    if(!(sym %in% rownames(profile))){
      prob[pos] = 0
    } else {
      prob[pos] = profile[sym,pos]
    }
  }
  #take the product of all positions
  kmer_prob = prod(prob)
  return(kmer_prob)
}

testthat::test_that("",{
  k = c("C","A", "B")
  q = cbind(c(0.25,0.5,0.25),c(0.5,0.5,0), c(0.25,0.5,0.25))
  rownames(q) = c("A","B","C")
  output = kmer_prob(k, q)
  
  ans = 0.25*0.5*0.5
  testthat::expect_equal(ans, output)
})

testthat::test_that("",{
  #case where kmer has letter not present in profile
  k = c("C","A", "D")
  q = cbind(c(0.25,0.5,0.25),c(0.5,0.5,0), c(0.25,0.5,0.25))
  rownames(q) = c("A","B","C")
  output = kmer_prob(k, q)
  ans = 0
  testthat::expect_equal(ans, output)
})

testthat::test_that("",{
  #case where kmer has letter not present in profile
  k = c("C","A", "B","B")
  q = cbind(c(0.25,0.5,0.25),c(0.5,0.5,0), c(0.25,0.5,0.25), c(0.2,0.1,0.7))
  rownames(q) = c("A","B","C")
  output = kmer_prob(k, q)
  ans = 0.25*0.5*0.5*0.1
  testthat::expect_equal(ans, output)
})
