#######################################################
#	We'll specify a piecewise constant function f(x) 
# with two vectors: a list of (n+1) boundaries 
#
#	b_1 < b_2 < ... < b_{n+1}
#
# and a list of values {v_1, v_2, ..., v_n}.
# The convention will be that if x lies in
# the j-th half-open interval [b_j, b_{j+1}),
# then f(x) = v_j. 
#
# We then want to write a function that
# accepts a new list of boundaries 
#
#	c_1 < ... c_{m+1}
#
# and computes an associated list of values 
# {u_1, ..., u_m} that are means of the original 
# step function over the new intervals. That is,
#
#	u_j = (1/(c_{j+1} - c_j)) *\int_{c_j}^c_{j+1} f(x) dx.
#####################################################

do.tests <- FALSE

##################################################################################
#	Integrate a step function over an interval
##################################################################################

#lb = lower bound (numeric vector)
#ub = upper bound (numeric vector)
#sf = step function, list of 2 numeric vectors (bdy,val)

integrate.step.func <- function( lb, ub, sf ){
	# If need be, swap the order of the bounds
	if( ub < lb ) {
		# Swap the bounds
		tmp = lb ; lb = ub ; ub = tmp ;
		
		# Note the change of sign
		my.sign = -1.0
	}
	else {
		my.sign = 1.0
	}
	
	# Deal with two trivial cases
	n.bins <- length( sf$val )
	if( (ub <= sf$bdy[1]) || (lb >= sf$bdy[n.bins+1]) ) {
		# The interval (ub,lb) doesn't overlap the support of sf.
		return( 0.0 )
	}
	
	# Find the indices of the bins containing ub and lb.
	j.lb <- bsearch.bins( lb, sf$bdy )
	if( j.lb == 0 ) {
		# Left part of the interval lies outside the support of sf
		lb <- sf$bdy[1]
		j.lb <- 1
	}
	
	j.ub <- bsearch.bins( ub, sf$bdy )
	if( j.lb == (n.bins + 1) ) {
		# Right part of the interval lies outside the support of sf
		ub <- sf$bdy[n.bins+1]
		j.ub <- n.bins
	}
	
	j <- j.lb
	integral <- 0.0
	while( j < j.ub ) {
		integral <- integral + (sf$bdy[j+1] - lb) * sf$val[j]
		lb <- sf$bdy[j+1]
		j <- j + 1
	}
	
	# Add the final contribution and return
	integral <- integral + (ub - lb) * sf$val[j.lb]
	return( my.sign * integral )
}

##################################################################################
#	Check whether a step function appears to be set up correctly
##################################################################################

# Are the entries in x in strictly increasing order?
is.increasing <- function( x ) {
	n <- length(x)
	return( all(x[1:(n-1)] < x[2:n]) )
}

# Does a step function appear to be sensibly defined?
validate.step.func <- function( sf ) 
{
	# Check that we have one value and 
	# two boundary points per interval
	n.intervals <- length(sf$val)
	if( n.intervals == 0 ) {
		warning( "validate.step.func(): number of intervals is zero." )
		return( FALSE )
	}
	else if( length(sf$bdy) != (n.intervals + 1) ) {
		warning( "validate.step.func(): lengths of bdys and vals don't match." )
		return( FALSE )
	}
	
	# Check that the boundary points are arranged
	# in strictly increasing order.
	if( !is.increasing( sf$bdy ) ) {
		warning( "validate.step.func(): list of bdy pts isn't strictly increasing." )
		return( FALSE )
	}
	
	# If we reach this line, all is well.
	return( TRUE )
}

if( do.tests ) {
	# Do testing: all the tests should come up TRUE
	verdicts <- c(
		is.increasing( c(1, 2) ),
		!is.increasing( c(1, 2, 2, 3) ),
		!is.increasing( c(1, 3, 2) ),
		validate.step.func( list(bdy=c(1,2,3), val=c(0.5, 1.5)) ),
		!validate.step.func( list(bdy=c(1,3,2), val=c(0.5, 1.5)) ),
		!validate.step.func( list(bdy=c(1,2), val=c(0.5, 1.5)) )
	)
	
	all( verdicts )
}

#some are tests are meant to go off for testing reasons

##################################################################################
# Given a list of N bin-boundaries and a value x, return the index of the bin
# that contains x. If x lies to the left of all the bins, return 0, and if
# x is to the right of all the bins, return (N+1). 
#
# N.B. The bins are assumed to be half-open and closed at the left.
##################################################################################

bsearch.bins <- function( x, bdy ) {
	# Deal with special cases first
	n.bins <- length(bdy) - 1
	if( x < bdy[1] ) {
		# x lies to the left of all bins
		return( 0 )
	}
	else if( x == bdy[1] ) {
		# x lies on the leftmost, closed boundary of the bins
		return( 1 )
	}
	else if( x >= bdy[n.bins + 1] ) {
		# x lies to the right of all bins, or on their rightmost (so, open) boundary.
		return( n.bins + 1 )
	}
	
	# If we reach these lines, bdy[1] < x < bdy[n.bins+1]
	left = 1
	right = n.bins + 1
	gap = n.bins
	while( gap > 1 ) {
		mid = left + (gap %/% 2) # integer division
		if( bdy[mid] < x ) {
			left = mid
		} else if( bdy[mid] == x ) {
			# x lies on a bin boundary, so gets assigned to a left endpoint
			return( mid )
		}
		else { # x < bdy[mid] 
			right = mid
		}
			
		gap = right - left 
	}
	
	# If we reach this line, bdy[left] < x < bdy[left+1]
	return( left )
}

if( do.tests ) {
	# Do testing: all the tests should come up TRUE
	bdys <- c( 1:5 )
	verdicts <- c(
		bsearch.bins( 0.0, bdys ) == 0,
		bsearch.bins( 1.0, bdys ) == 1,
		bsearch.bins( 1.3, bdys ) == 1,
		bsearch.bins( 2.0, bdys ) == 2,
		bsearch.bins( 4.5, bdys ) == 4,
		bsearch.bins( 5.0, bdys ) == 5,
		bsearch.bins( 5.5, bdys ) == 5
	)
	
	all( verdicts )
}

##################################################################################
# Finally, the main function
##################################################################################

#input: orig.sf : original step function. This is a list of 2 numeric vectors, bdy and val.
#new.bdy : A numeric vector containing the new boundary points. 

#output: New step function with resampled values. List of 2 numeric vectors, bdy and val. 

resample.step.func <- function( orig.sf, new.bdy ) {
	# Check the args
	if( !validate.step.func(orig.sf) ) {
		warning( "resample.step.func(): orig.sf is not valid." )
		return( NA )
	}
	else if( length(new.bdy) == 1 ) {
		warning( "resample.step.func(): list of new bdys pts has only one element." )
		return( FALSE )
	}
	else if( !is.increasing(new.bdy) ) {
		warning( "resample.step.func(): list of new bdys pts isn't strictly increasing." )
		return( FALSE )
	}
	
	# Collect some useful information
	n.orig <- length( orig.sf$val )
	n.new <- length( new.bdy ) - 1

	# Initialise the result
	new.sf <- list( bdy=new.bdy, val=rep( 0.0, n.new) )
	
	# Deal with two trivial cases where the sets of intervals are disjoint
	if( orig.sf$bdy[n.orig+1] <= new.bdy[1] ) {
		warning( "resample.step.func(): new and original intervals don't overlap." )
		return( new.sf )
	}
	else if( new.bdy[n.new+1] <= orig.sf$bdy[1] ) {
		warning( "resample.step.func(): new and original intervals don't overlap." )
		return( new.sf )
	}
	
	# The main event: integrate the original step function
	# across the intervals of the new one and then compute
	# the mean values
	for( j in 1:n.new ) {
		integral <- integrate.step.func( new.bdy[j], new.bdy[j+1], orig.sf )
		new.sf$val[j] <- integral / (new.bdy[j+1] - new.bdy[j])
	}
	
	return( new.sf )
}

if( do.tests ) {
	# Define a function to check whether two vectors agree
	vecs.agree <- function( u, v ) {
		if( length(u) != length(v) ) {
			verdict <- FALSE
		}
		else {
			abs.diff <- sum(abs(u - v))
			tol <- .Machine$double.eps * (sum(abs(u)) + sum(abs(v)))
			verdict <- abs.diff < tol
		}
		
		return( verdict )
	}
	
	orig.sf <- list( bdy=(0:10)/10, val=rep(1, 10))
	new.sf1 <- resample.step.func( orig.sf, (0:2)/2)
	vecs.agree( new.sf1$val, c(1,1))
	
	new.sf2 <- resample.step.func( new.sf1, (0:20)/20 )
	vecs.agree( new.sf2$val, rep(1,20) )
}


