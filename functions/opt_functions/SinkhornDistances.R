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
    # that p has only nonzero entries.
    engine.args <- list(
    	q=q, epsilon=epsilon, tol=tol, 
    	max.cycles=max.cycles, cycles.per.check=cycles.per.check,
    	print.progress=print.progress
    )
    
    p.has.zeros <- any(p == 0.0)
    if( p.has.zeros ) {
        idx <- p > 0
        engine.args$p <- p[idx]
        engine.args$cost.mat <- cost.mat[,idx]
    } else {
        engine.args$p <- p
        engine.args$cost.mat <- cost.mat
    }
    
    # Invoke either the scalar or sequence versions of the Sinkhorn algorithm.
    if( length(epsilon) == 1 ) {
    	result <- RcppSinkhornEngine( engine.args )
    } else {
    	# We have a sequence of epsilon values
    	result <- RcppSinkhornEpsSeqEngine( engine.args )
    }

    # If p contained zeroes, adjust the result accordingly
    if( p.has.zeros ) {
    	tmp <- matrix(rep( 0.0, length(p)*length(q) ), nrow=length(q) )
		for( j in 1:length(epsilon)) {
			tmp[,idx] <- result$transport.mat[[j]]
			result$transport.mat[[j]] <- tmp
		}
    }

    return( result )
}
