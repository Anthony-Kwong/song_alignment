####################################################################
#	Test the Sinkhorn engine on a simpler problem
####################################################################

library( RColorBrewer )
pair.pal <- brewer.pal( 12, "Paired" )

# Load all of our tools
do.tests <- TRUE
source( "SinkhornDistances.R" )
source( "ParabolicPdf.R" )

#################################################################################
# Make histograms of two different distribs, one the original parabolic
# distib and another that's a mixture of two parabolic bumps.
#################################################################################
n.bins <- 50
x.vals <- seq( from=-1, to=1, length.out=(n.bins + 1) )
dx = x.vals[2] - x.vals[1] 

p.mid.pts <- 0.5 * (x.vals[(1:n.bins) + 1] + x.vals[1:n.bins])
p.hist <- sapply(x.vals[(1:n.bins) + 1], parabolic.cdf) - sapply(x.vals[1:n.bins], parabolic.cdf)

q.mix.alpha <- 0.25
q.offset <- 2.2 
q.mid.pts <- c( p.mid.pts - q.offset, p.mid.pts + q.offset )
q.hist <- c( q.mix.alpha * p.hist, (1.0 - q.mix.alpha) * p.hist ) 

#################################################################################
# Plot everything: this is surprisingly tedious
#################################################################################

# Define a function to plot one of the bumps
add.bump <- function( x, y, fill.color, edge.color ) {
	# Plot a filled polygon
	poly.x <- rep( x, each=2 )
	poly.y <- c( 0.0, rep(y, each=2), 0.0)
	polygon( poly.x, poly.y, col=fill.color )
	
	# Add a border along the top
	lines( poly.x, poly.y, col=edge.color )
}

# Begin with empty axes
plot( x.vals, c(p.hist, 0), type="n", 
    xlim=c(-3,3), ylim=c(0, 0.8),
    xlab="x", ylab="Density", main="A toy optimal transport problem"
)

# Add the three bumps
add.bump( x.vals, p.hist/dx, pair.pal[1], pair.pal[2] )
add.bump( x.vals - q.offset, (q.mix.alpha/dx) * p.hist, pair.pal[3], pair.pal[4] )
add.bump( x.vals + q.offset, ((1.0 - q.mix.alpha)/dx) * p.hist, pair.pal[3], pair.pal[4] )

#################################################################################
# Build a cost matrix using an L_k norm
#################################################################################

cost.mat.k <- 2 # k=2 is the Wasserstein metric
mat.tmp <- rep( 0.0, length(p.mid.pts)*length(q.mid.pts) )
cost.mat <- matrix( mat.tmp, ncol=length(p.mid.pts) )
for( j in 1:n.bins ) {
	cost.mat[,j] <- abs(q.mid.pts - p.mid.pts[j])^cost.mat.k
}

################################################################
# Compute optimal transport plans and Sinkhorn distances
# for a range of values of epsilon.
################################################################

eps.vals <- 1.0/2^(1:7)

n.eps <- length( eps.vals )
dist.vals <- rep( 0.0, n.eps )
entropy.vals <- rep( 0.0, n.eps )
elapsed <- rep( 0.0, n.eps )
n.cycles <- rep( 0, n.eps )
verdicts <- rep( FALSE, n.eps )
for( j in 1:n.eps ) {
	sinkhorn.eps <- eps.vals[j]
	optimal.plan <- entropyRegularisedKOT( 
		p.hist, q.hist, cost.mat, sinkhorn.eps, 
		max.cycles=2000, tol=0.01*max(q.hist),
		print.progress=TRUE
	)
	
	verdicts[j] <- optimal.plan$converged
	n.cycles[j] <- optimal.plan$n.cycles
	dist.vals[j] <- (optimal.plan$optimal.cost)^(1.0/cost.mat.k)
	entropy.vals[j] <- optimal.plan$entropy
	elapsed[j] <- optimal.plan$elapsed
}

################################################################
# Plot the plan associated with the smallest epsilon
################################################################

# Add a block of zeroes into the middle of the optimal plan
# to account for those bins where nothing goes.
row.mid.pts <- c(
	p.mid.pts - q.offset,
	seq( from=(-q.offset + 1) + dx/2, to=(q.offset - 1) - dx/2, by=dx ),
	p.mid.pts + q.offset
)

n.extra.rows <- length(row.mid.pts) - 2*n.bins 
transport.plan <- rbind(
	optimal.plan$transport.mat[1:n.bins,],
	matrix( rep(0.0, n.extra.rows*n.bins), ncol=n.bins ),
	optimal.plan$transport.mat[(n.bins+1):(2*n.bins),]
)

image(
	y=p.mid.pts,
	x=row.mid.pts,
	z=transport.plan,
	ylab="Source distrib.",
	xlab="Target distrib."
)

# Print a summary
summary.df <- data.frame(
	Epsilon = eps.vals,
	Distance = dist.vals,
	Entropy = entropy.vals,
	Converged = verdicts,
	N.Cycles = n.cycles,
	Elapsed.Millisecs = elapsed
)
	
summary.df 
sum(elapsed)

