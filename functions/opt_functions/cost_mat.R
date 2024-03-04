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
cost_mat <- function(p,q, k = 2){
  #get vectors for source distribution
  
  #frequency
  f = p$freq
  #same sampling rate for every recording so these should be equal
  testthat::expect_equal(p$freq, q$freq)
  
  #time
  tp = p$time
  tq = q$time
  
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
  
  C = matrix(0, ncol = plen, nrow = qlen)
  #loop over rows
  for(i in 1:qlen){
    #loop over cols
    for(j in 1:plen){
      Q = q_vals[i,]
      P = p_vals[j,]
      C[i,j] = (Q$freq.ticks - P$freq.ticks)^k + a^(k)*(Q$time.ticks - P$time.ticks)^k
    }
  }
  
  return(C)
}

#test


# cost.mat.k <- 2 # k=2 is the Wasserstein metric
# mat.tmp <- rep( 0.0, length(p.mid.pts)*length(q.mid.pts) )
# cost.mat <- matrix( mat.tmp, ncol=length(p.mid.pts) )
# for( j in 1:n.bins ) {
#   cost.mat[,j] <- abs(q.mid.pts - p.mid.pts[j])^cost.mat.k
# }
