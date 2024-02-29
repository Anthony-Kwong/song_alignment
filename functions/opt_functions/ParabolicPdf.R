##############################################
#	Set up and check a parabolic PDF
##############################################

parabolic.pdf <- function( x ) { ifelse( abs(x) <= 1, 0.75 * (1 - x^2), 0 ) }
parabolic.cdf <- function( x ) { 
	if( x < -1 ) { return(0.0) }
	else if( x <= 1 ) { return( 0.75 * (x - x^3/3.0 + 2.0/3.0) ) }
	else { return( 1.0 ) }
}

invert.parabolic.cdf <- function( u )
{
	if( (u < 0) || (u > 1) ) { return(NA) }
	else if( u == 0 ) { return( -1 ) }
	else if( u == 1 ) { return( 1 )  }
	else {
		f <- function( x ) { parabolic.cdf(x) - u }
		soln <- uniroot( f, lower=-1, upper= 1)
		return( soln$root )
	}
}

random.parabolic <- function( n, mean=0, w=1 ) 
{
	x <- sapply( runif(n), invert.parabolic.cdf )
	return( w*x + mean )
}

if( do.tests ) {
	x.vals <- seq( from=-1.5, to = 1.5, length.out=301 ) ;
	y.vals <- sapply( x.vals, parabolic.pdf )
	plot( x.vals, y.vals, type="l", xlab="x", ylab="Density")

	hist( random.parabolic(10000), breaks="FD", xlab="x", ylab="Density", freq=FALSE )
	lines( x.vals, y.vals )
	
	integrate( parabolic.pdf, lower=-1, upper=1 ) 
}
