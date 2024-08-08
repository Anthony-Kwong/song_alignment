library( alabama )

do.tests <- TRUE
check.gradients <- FALSE 

###########################################################################
#	Given father-to-son transmission pairs (x_i, y_i), find a matrix T 
# that minimises a certain natural objective function.
###########################################################################

estimateTransmissionMatrix <- function( father.mat, son.mat, cost.mat=NULL, tol=1.0e-10 ) {
	# Check the args
	if( (nrow(father.mat) != nrow(son.mat)) || (ncol(father.mat) != ncol(son.mat)) ) {
		warning( "Father and son matrices have different shapes." )
		return( Null )
	} else if( any(father.mat < 0.0) || any(father.mat > 1.0) ) {
		warning( "Father matrix has entries outside [0.0, 1.0]" )
		return( Null )
	} else if( any(son.mat < 0.0) || any(son.mat > 1.0) ) {
		warning( "Son matrix has entries outside [0.0, 1.0]" )
		return( Null )
	} else if( max(abs(colSums(father.mat) - 1.0)) > tol ) {
		warning( "Columns of father matrix don't sum to 1.0'" )
		return( Null )
	} else if( max(abs(colSums(son.mat) - 1.0)) > tol ) {
		warning( "Columns of son matrix don't sum to 1.0'" )
		return( Null )
	}
	
	# Extract some useful sizes
	n.notes <- nrow(son.mat)
	n.pairs <- ncol(son.mat)
	T.len <- n.notes*n.notes
	
	# Build the objective function and its gradient
	objective <- function( T.vec ) {
		T.mat <- matrix( T.vec, nrow=n.notes )
		predicted.son.mat <- T.mat %*% father.mat
		diff.mat <- son.mat - predicted.son.mat
		return( 0.5 * sum(diff.mat * diff.mat) )
	}
	
	grad.objective <- function( T.vec ) {
		T.mat <- matrix( T.vec, nrow=n.notes )
		predicted.son.mat <- T.mat %*% father.mat
		grad.as.mat <- matrix( rep( 0.0, T.len ), nrow=n.notes )
		for( i in 1:n.pairs ) {
			predicted.son.vec <- T.mat %*% father.mat[,i]
			diff.vec <- as.vector(predicted.son.vec - son.mat[,i])
			grad.as.mat <- grad.as.mat + outer(diff.vec, father.mat[,i])
		}
		
		return( as.vector(grad.as.mat) )
	}

	# Build the inequality constraints and their gradients (a Jacobian). 
	# The constraints need to be strict, so we'll demand 0 < T_{i,j} < 1.
	inequality.constraints <- function( T.vec ) {
		h <- c( T.vec, 1.0 - T.vec )
		return( h )
	}
	
	# The Jacobian is constant, so we compute it once and for all.
	inequality.jac.mat <- rbind( diag(1.0, nrow=T.len), diag(-1.0, nrow=T.len))
	inequality.jacobian <- function( T.vec ) {
		return( inequality.jac.mat  )
	}
	
	# Build the equality constraints and their gradients (a Jacobian). 
	# We require the columns of T to sum to 1.0.
	equality.constraints <- function( T.vec ) {
		T.mat <- matrix( T.vec, nrow=n.notes )
		h <- colSums( T.mat ) - 1.0 
		return( h )
	}
	
	# The Jacobian is constant, so we compute it once and for all.
	for( j in 1:n.notes ) {
		tmp <- rep( 0.0, n.notes )
		tmp[j] <- 1.0 
		crnt.row <- rep( tmp, each=n.notes )
		if( j == 1 ) {
			equality.jac.mat <- crnt.row
		}
		else {
			equality.jac.mat <- rbind( equality.jac.mat, crnt.row )
		}
	}
	
	equality.jacobian <- function( T.vec ) {
		return( equality.jac.mat )
	}
	
	# Choose an initial state that favours preservation of the father's note.
	delta <- 0.1 # fraction of notes that change
	off.diag.entry <- delta / (n.notes - 1)
	T.mat.zero <- matrix( rep(off.diag.entry, T.len), nrow=n.notes )
	T.mat.zero <- T.mat.zero + diag( (1.0 - delta - off.diag.entry), nrow=n.notes )
	T.vec.zero <- as.vector( T.mat.zero )
	
	# If we're testing, check all the gradients
	if( do.tests && check.gradients ) {
		# Test all the gradients
		delta.t <- 1.0e-4
		
		# First the objective
		computed.grad <- grad.objective( T.vec.zero )
		fd.grad <- rep( 0.0, T.len )
		for( j in 1:T.len ) {
			T.plus <- T.vec.zero
			T.plus[j] <- T.plus[j] + delta.t
			obj.plus <- objective( T.plus )
			
			T.minus <- T.vec.zero
			T.minus[j] <- T.minus[j] - delta.t
			obj.minus <- objective( T.minus )
			
			fd.grad[j] <- (obj.plus - obj.minus) / (2.0 * delta.t)
		}
		
		max.abs.diff <- max(abs(computed.grad - fd.grad))
		print( sprintf("grad.objective: %g",  max.abs.diff), quote=FALSE )
		
		# Then the inequality constraints
		computed.jac <- inequality.jacobian( T.vec.zero )
		for( i in 1:nrow(computed.jac) ) {
			computed.grad <- computed.jac[i,]
			fd.grad <- rep( 0.0, T.len )
			for( j in 1:T.len ) {
				T.plus <- T.vec.zero
				T.plus[j] <- T.plus[j] + delta.t
				hin.plus <- inequality.constraints( T.plus )[i]
				
				T.minus <- T.vec.zero
				T.minus[j] <- T.minus[j] - delta.t
				hin.minus <- inequality.constraints( T.minus )[i]
				
				fd.grad[j] <- (hin.plus - hin.minus) / (2.0 * delta.t)
			}
			
			max.abs.diff <- max(abs(computed.grad - fd.grad))
			print( sprintf("inequality.jacobian[%d,]: %g", i, max.abs.diff), quote=FALSE )
		}
		
		# And finally, the equality constraints
		computed.jac <- equality.jacobian( T.vec.zero )
		for( i in 1:nrow(computed.jac) ) {
			computed.grad <- computed.jac[i,]
			fd.grad <- rep( 0.0, T.len )
			for( j in 1:T.len ) {
				T.plus <- T.vec.zero
				T.plus[j] <- T.plus[j] + delta.t
				heq.plus <- equality.constraints( T.plus )[i]
				
				T.minus <- T.vec.zero
				T.minus[j] <- T.minus[j] - delta.t
				heq.minus <- equality.constraints( T.minus )[i]
				
				fd.grad[j] <- (heq.plus - heq.minus) / (2.0 * delta.t)
			}
			
			max.abs.diff <- max(abs(computed.grad - fd.grad))
			print( sprintf("equality.jacobian[%d,]: %g", i, max.abs.diff), quote=FALSE )
		}
	}

	# Do the optimisation and return
	result <- constrOptim.nl(
		T.vec.zero,
		fn=objective, gr=grad.objective,
		hin=inequality.constraints, hin.jac=inequality.jacobian,
		heq=equality.constraints, heq.jac=equality.jacobian,
		control.outer=list( trace=FALSE, mu0=1.0e-9, sig0=1.0e-9 ), 
		control.optim=list( abstol=1.0e-12 )
	)
	
	if( result$convergence != 0 ) {
		print( result$message )
	}
	
	return( result )
}

if( do.tests ) {
	library( gtools )  # for rdirichlet()
	
	# Set the sizes
	n.notes <- 12
	n.pairs <- 1000
	T.len <- n.notes * n.notes
	
	# Generate synthetic data and test the optimiser
	alpha.vec <- rep( 1.0, n.notes ) # For a dirichlet distribution
	T.mat <- t(rdirichlet( n.notes, alpha.vec))
	
	alpha.vec <- rep( 0.25, n.notes ) # Favours small entries
	father.mat <- t(rdirichlet( n.pairs, alpha.vec))
	son.mat <- T.mat %*% father.mat
	
	# Add a bit of noise to son.mat. This is a bit of a faff as
	# we need the columns of the result to be non-negative and sum to 1.0
	eps <- 0.005
	noise.vals <- runif( n.pairs*n.notes, min=-eps, max=eps)
	noise.mat <- matrix( noise.vals, nrow=n.notes )
	son.mat <- son.mat + noise.mat
	son.mat <- ifelse( son.mat < 0, 0, son.mat )
	son.mat <- ifelse( son.mat > 1, 1, son.mat )
	son.mat <- son.mat %*% diag( 1.0/colSums(son.mat) )
	
	# Do the thing
	result <- estimateTransmissionMatrix( father.mat, son.mat )
	fitted.T.mat <- matrix( result$par, nrow=n.notes )
	colSums( fitted.T.mat )
	T.mat - fitted.T.mat
}

#use on real data

sf.mat = readRDS(file ="./results/seq_evo/sf_mat.rds")
son.mat = readRDS(file ="./results/seq_evo/son_mat.rds")
res = estimateTransmissionMatrix(father.mat = sf.mat, son.mat = son.mat)

matrix(res$par, ncol = 16)
