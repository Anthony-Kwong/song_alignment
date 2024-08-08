#' note_compare function
#' 
#' Computes sinkorn distances between a pair of note spectrograms.
#'
#' @param n1 : Spectrogram of note1 as a list.
#' @param n2 : Spectrogram of note2 as a list.
#' @param a: Numeric. Weight to modify the time dimension. 
#' @param epsilon : Numeric. Epsilon value to use for optimal transport problem. (fast method)
#' @param max.cycles : Numeric. Max number of iterations to go through before exiting.
#' @param k : Numeric. Power for the Wasserstein distance.
#'
#' @return The output of entropyRegularisedKOT() from the SinkhornDistances.R.
#' @export
#'
#' @examples

#source dependencies
source("./SinkhornDistances.R")
source("./cost_mat.R")
#set tests
doTests = FALSE

note_compare <- function(n1,n2,a, epsilon=0.1, max.cycles=1000, k=2){
  
  #get amplitude matrices
  A1 = n1$amp
  A2 = n2$amp
  
  #turn amplitudes positive and normalise. 
  As = lapply(list(A1,A2), function(X){
    #turn everything positive, amplitudes are recorded from a baseline
    pos = X - min(X)
    #normalise
    pos/sum(pos)
  })
  
  #vectorise matrices for optimal transport ----
  vA1 = as.vector(As[[1]])
  vA2 = as.vector(As[[2]])
  
  #check correct values in vectorization
  if(doTests){
    D1 = As[[1]]
    testthat::expect_equal(D1[,1], vA1[1:nrow(D1)])
    testthat::expect_equal(D1[,2], vA1[(nrow(D1)+1):(2*nrow(D1))])
    D2 = As[[2]]
    testthat::expect_equal(D2[,1], vA2[1:nrow(D2)])
    testthat::expect_equal(D2[,3], vA2[(2*nrow(D2)+1):(3*nrow(D2))])
  }
  
  #Cost matrix
  C = cost_mat(p = n1, q = n2, a = a, k)
  if( is.null(C) ) {
    return( NULL )
  }
  
  #compute sinkorn distances
  res = entropyRegularisedKOT(p = vA1, q = vA2, cost.mat = C, epsilon=epsilon, max.cycles=max.cycles )
  
  return(res)
}
