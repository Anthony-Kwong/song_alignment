#' margin_opt function
#' 
#' Marginalise the distributions of 2 notes and computes Sinkorn distances. We marginalise over
#' time and frequency.
#'
#' @param n1 : Spectrogram of note1 as a list.
#' @param n2 : Spectrogram of note2 as a list.
#'
#' @return : A list of the output of entropyRegularisedKOT() from the SinkhornDistances.R. One for 
#' time dimension, one for frequency dimension.
#' @export
#'
#' @examples
library(Rcpp)
sourceCpp("./functions/opt_functions/costmat_vec.cpp")
source("./functions/opt_functions/SinkhornDistances.R")
margin_opt <- function(n1,n2){
  #take marginals
  
  #processing amplitude matrices
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
  
  #frequency ----
  
  #frequency marginal 
  f_n1 = rowSums(As[[1]])
  f_n2 = rowSums(As[[2]])
  
  #get dimensions
  freq_dim = length(n1$freq)
  if(freq_dim != length(n2$freq)){
    stop("frequency dimensions not matching")
  }
  
  #get frequency vector
  if(setequal(n1$freq,n2$freq) ==F ){
    stop("frequency vectors not matching")
  }
  freq = n1$freq
  
  #cost matrix
  C = costmat_vec(freq,freq)
  
  #compute OPT
  f_opt = entropyRegularisedKOT(p = f_n1, q = f_n2, cost.mat = C)
  
  #time ----
  At_1 = colSums(As[[1]])
  At_2 = colSums(As[[2]])
  
  #time: make each interval symmetric about zero
  t_n1 = n1$time
  t_n2 = n2$time
  tp = t_n1 - (t_n1[1] + t_n1[length(t_n1)]) / 2
  tq = t_n2- (t_n2[1] + t_n2[length(t_n2)]) / 2
  #cost matrix
  C = costmat_vec(tp,tq)
  #compute OPT
  t_opt = entropyRegularisedKOT(p = At_1, q = At_2, cost.mat = t(C))
  
  output = list(freq = f_opt, time = t_opt)
  return(output)
}

#add tests
