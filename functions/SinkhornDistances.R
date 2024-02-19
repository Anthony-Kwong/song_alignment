entropyRegularisedKOT <- function( p, q, cost.mat, gamma=1.0, tol=1.0e-8, max.cycles=500, tiny.val=1.0e-200 ) {
    # Check the args
    if( gamma <= 0.0 ) {
        warning( "entropyRegularisedKOT: gamma should be positive" ) 
        return( NA )
    } else if( any( p < 0 ) ) {
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

    # Drop any entries in p that are zero and drop the corresponding columns
    # of the cost matrix.
    p.has.zeros <- any(p == 0.0)
    if( p.has.zeros ) {
        idx <- p > 0
        nonzero.p <- p[idx] 
        cost.mat <- cost.mat[,idx]
    } else {
        nonzero.p <- p
    }
   
    # The iterative algorithm below alternately adjusts the digonal matrices D_1 and D_2 
    # whose digonal entries are stored in vectors u an v, respectively.
    n.cycles <- 0
    exp.cost.mat <- exp( -cost.mat / gamma ) # elementwise exponentiation
    v <- rep( 1, length(nonzero.p) ) # Initial value for diagonal entries of D_2
    KOT.mat <- q %o% p # An outer product that gives the solution when gamma tends to infinity
    while( (n.cycles == 0) || ((max.diff > tol) && (n.cycles < max.cycles)) ) {
        # Update the two diagonal matrices, taking care to avoid zeroes
        u <- q / exp.cost.mat %*% v
        u <- ifelse( u < tiny.val, tiny.val, u ) 
        
        v <- nonzero.p / t(exp.cost.mat) %*% u
        v <- ifelse( v < tiny.val, tiny.val, v )  

        # Recompute the solution and the convergence score
        uKv <- diag( as.vector(u) ) %*% exp.cost.mat %*% diag( as.vector(v) )
        max.diff <- max(abs(KOT.mat - uKv)) 
        
        # On to the next cycle
        KOT.mat <- uKv
        n.cycles <- n.cycles + 1 
    }

    # If p contains zeroes, assemble a KOT matrix with 
    # columns of zeroes
    if( p.has.zeros ) {
        tmp <- matrix(rep( 0.0, length(p)*length(q) ), nrow=length(q) )
        tmp[,idx] <- KOT.mat
        KOT.mat <- tmp
    }
    
    sinkhorn.dist <- sum( KOT.mat * cost.mat ) # Product here is elementwise
    return( list( transport.mat=KOT.mat, sinkhorn.dist=sinkhorn.dist, gamma=gamma, converged=(max.diff < tol), n.cycles=n.cycles ) )
}