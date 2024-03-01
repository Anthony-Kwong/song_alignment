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
  
}

# cost.mat.k <- 2 # k=2 is the Wasserstein metric
# mat.tmp <- rep( 0.0, length(p.mid.pts)*length(q.mid.pts) )
# cost.mat <- matrix( mat.tmp, ncol=length(p.mid.pts) )
# for( j in 1:n.bins ) {
#   cost.mat[,j] <- abs(q.mid.pts - p.mid.pts[j])^cost.mat.k
# }
