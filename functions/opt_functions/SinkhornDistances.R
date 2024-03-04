####################################################################
#	Here is a function that computes entropy-regularised 
# optimal-transport distances in the style of
#
#	Mario Cuturi (2013), Sinkhorn Distances: Lightspeed Computation 
#	of Optimal Transportation Distances, Advances in Neural Information 
#	Processing Systems 26:2292–2300.
#
# mrm: ManUni & Whalley Range, 13-16 Feb 2024.
####################################################################
library(Rcpp)
library(RcppArmadillo)

# Compile the C++ engine
sourceCpp( "SinkhornDistanceEngine.cpp" )

#for more theory see: https://lucyliu-ucsb.github.io/posts/Sinkhorn-algorithm/

####################################################################
# Checks the args and packages data for the engines.
####################################################################
#' Compute Sinkorn distances between 2 distributions.
#'
#' @param p : Numeric vector representing source distribution. 
#' @param q : Numeric vector representing destination distribution.
#' @param cost.mat : Numeric matrix of costs.
#' @param epsilon : Numeric, coefficient of the regression term
#' @param tol : numeric, tolerance level for convergence
#' @param max.cycles : numeric, max number of iterations to cycle through
#' @param cycles.per.check : numeric, number of cycles to go through before checking transport plan (for efficiency)
#' @param print.progress : logical
#'
#' @return An object containing the optimal transport plan and cost. 
#' @export
#'
#' @examples p = c(0.25, 0.25, 0.4, 0.1)
#' q = c(0.1, 0.15, 0.5, 0.1, 0.15)
#' C = matrix(1, nrow = 5, ncol = 4)
#' entropyRegularisedKOT(p=p,q=q,cost.mat = C)
#' 

####################################################################
# Checks the args and packages data for the engines.
####################################################################
entropyRegularisedKOT <- function( 
	p, q, cost.mat, epsilon=1.0, 
	tol=1.0e-6, max.cycles=200, cycles.per.check=10,
	print.progress=FALSE )
{
    # Check that epsilon is positive
    if( any(epsilon <= 0) ) {
        warning( "entropyRegularisedKOT: epsilon should be positive" ) 
        return( NA )
    } 
    
    # Check that p and q look like probability distributions
    if( any( p < 0 ) ) {
        warning( "entropyRegularisedKOT: p contains negative entries" ) 
        return( NA )
    } else if( abs(sum(p) - 1.0) > tol ) {
        warning( "entropyRegularisedKOT: the entries in p do not sum to one" ) 
        return( NA )
    } else if( any( q < 0 ) ) {
        warning( "entropyRegularisedKOT: q contains negative entries" ) 
        return( NA )
    } else if( abs(sum(q) - 1.0) > tol ) {
        warning( "entropyRegularisedKOT: the entries in q do not sum to one" ) 
        return( NA )
    }
		
    # Assemble the input data for the engines. The main work here is to ensure
    # that p and q have only nonzero entries.
    engine.args <- list(
    	epsilon=epsilon, tol=tol, 
    	max.cycles=max.cycles, cycles.per.check=cycles.per.check,
    	print.progress=print.progress
    )
    
	p.idx <- p > 0.0
	q.idx <- q > 0.0    
    engine.args$p <- p[p.idx]
    engine.args$q <- q[q.idx]
    engine.args$cost.mat <- cost.mat[q.idx, p.idx]
    
    # Invoke either the scalar or sequence versions of the Sinkhorn algorithm.
    if( length(epsilon) == 1 ) {
    	result <- RcppSinkhornEngine( engine.args )
    } else {
    	# We have a sequence of epsilon values
    	result <- RcppSinkhornEpsSeqEngine( engine.args )
    }

    # If p or q contained zeroes, adjust the result accordingly
	tmp <- matrix(rep( 0.0, length(p)*length(q) ), nrow=length(q) )
	if( length(epsilon) == 1 ) {
		tmp[q.idx,p.idx] <- result$transport.mat
		result$transport.mat <- tmp
	}
	else {
		for( j in 1:length(epsilon)) {
			tmp[q.idx,p.idx] <- result$transport.mat[[j]]
			full.results[[j]] <- tmp
		}
	}
	
    return( result )
}
