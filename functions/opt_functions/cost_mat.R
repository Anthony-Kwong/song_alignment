#' cost_mat function
#' 
#' Generates a cost matrix for the OPT problem of comparing notes. 
#'
#' @param p : List for the spectrogram of note1 (source distribution). Contains numeric vectors time,freq.
#' @param q : List for the spectrogram of note2 (destination distribution). Contains numeric vectors time,freq.
#' @param k : Numeric. Power of the distance function. Default k=2 for the Wasserstein metric.
#' @param a : Constant term coefficent. Default a=1.
#'
#' @return Numeric cost matrix. c = (f1-f2)^k + a*(t1-t2)^k. 
#' @export
#'
#' @examples
Rcpp::sourceCpp("./functions/opt_functions/costmat_C.cpp")
cost_mat <- function(p,q, a=1.0 ,k=2 ){
  #get vectors for source distribution
  
  #frequency
  f = p$freq
  #same sampling rate for every recording so these should be equal
  testthat::expect_equal(p$freq, q$freq)
  
  #time: make each interval symmetric about zero
  tp = p$time - (p$time[1] + p$time[length(p$time)]) / 2
  tq = q$time - (q$time[1] + q$time[length(q$time)]) / 2
  
  #get f and t values as tables ----
  vals = lapply(list(tp,tq), function(t){
    #index vectors to keep track of which elements to use
    p_reps = length(t)
    pf_index = rep(length(f):1, p_reps)
    pt_index = rep(seq(p_reps), each = length(f))
    #use indices to fetch the frequencies and times to use for computing costs
    pf = sapply(pf_index, function(x){f[x]})
    pt = sapply(pt_index, function(y){t[y]})
    tibble::tibble(freq.ticks = pf,time.ticks =pt)
  })
  
  p_vals = vals[[1]]
  q_vals = vals[[2]]
  
  #check correct dimensions
  testthat::expect_equal(nrow(p_vals), length(f)*length(p$time))
  testthat::expect_equal(nrow(q_vals), length(f)*length(q$time))
  
  #Compute matrix ----
  plen = nrow(p_vals)
  qlen = nrow(q_vals)
  
  C = costmat_C(pf = p_vals$freq.ticks, pt = p_vals$time.ticks, qf = q_vals$freq.ticks, qt = q_vals$time.ticks, a = a, k = k)
  
  return(C)
}

#add test to check correct inputs for costmat_C
f = seq(from = 0, to = 22, by = 4)
p = list(freq = f, time = seq(from=0,to=5, by = 3))
q = list(freq = f, time = seq(from=0,to=10, by = 4))

cost_mat()

