//#####################################################################
// 	Here is a C++ implementation of the a routine to calculate 
// entropy-regularised optimal-transport distances in the style of
//
//	Mario Cuturi (2013), Sinkhorn Distances: Lightspeed Computation 
//	of Optimal Transportation Distances, Advances in Neural Information 
//	Processing Systems 26:2292–2300.
//
// I wrote this as a replacement for an R implementation: both are 
// derived from a MATLAB version due to Cuturi and available at 
//
//	http://marcocuturi.net/SI.html
//
// as well as Remark 4.23 (Log-domain Sinkhorn) of
//
//	Gabriel Peyré and Marco Cuturi (2019), Computational Optimal Transport,
//	Foundations and Trends in Machine Learning 11(4-5):355-607.
//	https://arxiv.org/abs/1803.00567
// 
// mrm: Whalley Range, 15 Feb 2024
//#####################################################################
// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
# include <Rcpp/Benchmark/Timer.h>

# include <sstream> // string streams

using namespace Rcpp ;

//#####################################################################
//	The structure defined here holds workspaces and intermediate
// results. Its main use is to allow efficient re-use of results
// computed with a given value of epsilon as initial conditions for
// a computation with a smaller value of epsilon.
//#####################################################################

class SinkhornWorkspace 
{
	public:
		// These variables hold results that we'll want to copy
		arma::vec		f ;
		arma::vec		g ;
		arma::mat		P ;
		double			max_diff ;
		std::size_t		n_cycles ;
		
		// Construct an object to hold intermediate results used
		// in the computation of Sinkhorn distances.
		SinkhornWorkspace( 
			const arma::vec		&my_a,
			const arma::vec		&my_b,
			const arma::mat		&my_C,
			std::size_t			my_max_cycles,
			std::size_t			my_cycles_per_check,
			double				my_epsilon,
			double				my_tol
		)
		{
			// Copy the data that define the problem
			a = my_a ;	// target distrib
			b = my_b ;	// source distrib
			C = my_C ;	// cost matrix
			
			// Copy convergence data 
			tol = my_tol ;
			epsilon = my_epsilon ;
			max_cycles = my_max_cycles ;
			cycles_per_check = my_cycles_per_check ;
			converged = false ;	// for now
			
			// Precompute useful quantities
			eps_log_a = epsilon * arma::log(a) ;
			eps_log_b = epsilon * arma::log(b) ;
			neg_C_over_eps = (-1.0 / epsilon) * C ;
			
			// Preallocate storage for results
			f = arma::zeros( a.n_elem ) ; 
			g = arma::zeros( b.n_elem ) ;
			P.set_size(size(C)) ; 
			log_P.set_size(size(C)) ;
		
			// Preallocate various bits of scratch space
			Sfg.set_size(size(C)) ; 
			col_tmp.set_size( Sfg.n_rows ) ; 
			row_tmp.set_size( Sfg.n_cols ) ;
		}
		
		// The main thing we'll want to do is reset epsilon
		// and, at the same time, reset the convergence test.
		double get_epsilon( void ) const { return(epsilon) ; }
		void update_epsilon( double eps ) {
			// Be paranoid
			assert( epsilon > 0 ) ;
			assert( eps > 0.0 ) ;
			
			// Recompute epsilon-dependent quantities
			double old_eps = epsilon ;
			epsilon = eps ;
			
			f *= (epsilon/old_eps) ;
			g *= (epsilon/old_eps) ;

			eps_log_a *= (epsilon/old_eps) ;
			eps_log_b *= (epsilon/old_eps) ;
			neg_C_over_eps = (-1.0 / epsilon) * C ;
			
			// Whatever the previous verdict, we're starting afresh.
			converged = false ;
		}
		
		void get_P( arma::vec &my_P ) const { my_P = P ; }
		bool convergedQ( void ) const { return(converged) ; }
		
		// Entropy & cost for finished computations
		double optimal_cost( void ) const
			{ return( arma::accu( C % P ) ) ;  } // % is elementwise product 
			
		double entropy( void ) const // Eqn. (4.1) 
			{ return( -1.0 * arma::accu( (P % log_P) - P ) ) ;  }
			
		// The body of this member function, which is where the
		// main action occurs, appears below.
		void iterate( bool printProgress = false ) ;
			
	private:
		bool					converged ;
		double					tol ;
		double					epsilon ;
		std::size_t				max_cycles ;
		std::size_t				cycles_per_check ;
		arma::vec				a ; // Target distrib: see Eqn. (2.10)
		arma::vec				b ; // Source distrib: see Eqn. (2.10)
		arma::vec				eps_log_a ;
		arma::vec				eps_log_b ;
		arma::vec				col_tmp ;
		arma::vec				row_tmp ;
		arma::mat				C ;
		arma::mat				neg_C_over_eps ;
		arma::mat				Sfg ;
		arma::mat				log_P ;
} ;

//#####################################################################
// Prototypes for routines internal to this file.
//#####################################################################

static double soft_min( double epsilon, const arma::vec &z ) ;
static arma::vec &soft_mins_across_rows( arma::vec &result, double epsilon, const arma::mat &A ) ;
static arma::vec &soft_mins_down_cols( arma::vec &result, double epsilon, const arma::mat &A ) ;

//#####################################################################
// 	Given a vector (v_1, v_2, ..., v_n), compute
// 
//	-epsilon * log[ \sum_j exp(-z_j / epsilon) ].
//
// This quantity converges to the minimum of the z_j in the
// limit epsilon -> 0. See Remark 4.22 in Peyré and Cuturi.
//#####################################################################

// [[Rcpp::export]]
double soft_min( double epsilon, const arma::vec& z )
{
	// N.B. If some of the entries in z are zero, 
	// this function may not do what one expects.
	double sum_exp = arma::sum( arma::exp( -z / epsilon ) ) ;
	return( -epsilon * log(sum_exp) ) ;
}

//#####################################################################
// Two utilities defined at the bottom of page 78 in Peyré and Cuturi.
// N.B. These functions don't check their args.
//#####################################################################

arma::vec &soft_mins_down_cols( 
	arma::vec		&result,
	double			epsilon,
	const arma::mat	&A
) {
	for( std::size_t j=0 ; j < A.n_cols ; j++ ) {
		result(j) = soft_min( epsilon, A.col(j) ) ;
	}
	
	return( result ) ;
}

arma::vec &soft_mins_across_rows( 
	arma::vec		&result,
	double			epsilon,
	const arma::mat	&A
) {
	return( soft_mins_down_cols( result, epsilon, A.t() ) ) ;
}

//#####################################################################
// Do the heavy lifting involved in computing Sinkhorn distances. 
// N.B. This function does no error-checking: all that should be
// handled by some R-based wrapper. Also, p should have strictly
// nonzero entries.
//#####################################################################

// [[Rcpp::export]]
List RcppSinkhornEngine( const List &args )
{
	// Unpack the args and change notation to that used
	// by Peyré and Cuturi. At the same time, we construct
	// a workspace object that includes various useful
	// bits of pre-allocated storage.
	SinkhornWorkspace workspace( 
		as<arma::vec>(args["q"]), 
		as<arma::vec>(args["p"]),
		as<arma::mat>(args["cost.mat"]),
		args["max.cycles"], args["cycles.per.check"], 
		args["epsilon"], args["tol"] ) ;
		
	// Invoke the workhorse function and
	// do timing after the fashion of 
	// https://gallery.rcpp.org/articles/using-the-rcpp-timer/
	Timer timer ;
    timer.step("start") ;
	workspace.iterate( args["print.progress"] ) ;
	timer.step("stop") ;
	
	NumericVector res(timer) ;
	
    // Assemble the result and return
    List result ;
	result["transport.mat"] = workspace.P ;
	result["optimal.cost"] = workspace.optimal_cost() ;
	result["entropy"] = workspace.entropy() ;
	result["epsilon"] = args["epsilon"] ;
	result["elapsed"] = (res(1) - res(0)) / 1.0e6 ; // Rcpp's Timer does nanosecs
	result["converged"] = workspace.convergedQ() ;
	result["n.cycles"] = workspace.n_cycles ;
	result["log.row.scaling"] = workspace.f ;
	result["log.col.scaling"] = workspace.g ;
	
    return( result ) ;
}

//#####################################################################
// Handle a decreasing sequence of values of epsilon, using the result
// from the previous value as the starting point for the next.
//#####################################################################

// [[Rcpp::export]]
List RcppSinkhornEpsSeqEngine( const List &args )
{
	// Get the values of epsilon and ensure that 
	// they're in decreasing order.
	arma::vec epsilon = args["epsilon"] ;
	epsilon = arma::sort( epsilon, "descend" ) ;
	std::size_t	 n_eps = epsilon.n_elem ;
	
	// Pre-allocate the results, using R types
	List 			transport_mats(n_eps) ;
	List 			my_f(n_eps) ;
	List 			my_g(n_eps) ;
	NumericVector	cost(n_eps) ;
	NumericVector	entropy(n_eps) ;
	IntegerVector	n_cycles(n_eps) ;
	LogicalVector	converged(n_eps) ;
		
	// Unpack the args and change notation to that used
	// by Peyré and Cuturi. At the same time, we construct
	// a workspace object that includes various useful
	// bits of pre-allocated storage.
	SinkhornWorkspace workspace( 
		as<arma::vec>(args["q"]), 
		as<arma::vec>(args["p"]),
		as<arma::mat>(args["cost.mat"]),
		args["max.cycles"], args["cycles.per.check"], 
		epsilon[0], args["tol"] ) ;
	
	// Start a timer. 
	// We use a stringstream to set its labels
    Timer timer ;
    std::stringstream label ;
    
    timer.step( "start" ); 
	for( std::size_t j=0 ; j < n_eps ; j++ ) {
		// Update the value of epsilon
		workspace.update_epsilon( epsilon[j] ) ;
		
		// Invoke the workhorse function
		 workspace.iterate( args["print.progress"] ) ;
		 label << "epsilon = " << epsilon[j] ;
		 timer.step( label.str() ) ;
		 
		 // Copy the results into place
		 transport_mats[j] = workspace.P ;
		 cost[j] = workspace.optimal_cost() ;
		 entropy[j] = workspace.entropy() ;
		 n_cycles[j] = workspace.n_cycles ;
		 converged[j] = workspace.convergedQ() ;
		 my_f[j] = workspace.f ;
		 my_g[j] = workspace.g ;
	}
	
	// Get the timing results
	NumericVector res(timer) ;
	arma::vec elapsed(n_eps) ;
	for( std::size_t j=0 ; j < n_eps ; j++ ) {
		elapsed[j] = (res[j+1] - res[0]) / 1.0e6 ; // converts to millisecs
	}
	    
    // Assemble the result and return
    List result ;
	result["transport.mat"] = transport_mats ;
	result["optimal.cost"] = cost ;
	result["entropy"] = entropy ;
	result["epsilon"] = args["epsilon"] ;
	result["elapsed"] = elapsed ; 
	result["converged"] = converged ;
	result["n.cycles"] = n_cycles ;
	result["log.row.scaling"] = my_f ;
	result["log.col.scaling"] = my_g ;
	
    return( result ) ;
}

//#####################################################################
//	This is the workhorse: it manages the iterations of the 
// Sinkhorn algorithm.
//#####################################################################

void SinkhornWorkspace::iterate( bool printProgress ) 
//############################################################################
// Sinkhorn's algorithm, below, alternately updates two diagonal matrices,
// the logs of whose nonzero entries are stored in vectors f and g.
// Following Peyré and Cuturi's Section 4.4, Stability and Log-Domain Computations,
// we work with logs throughout, though the updates turn out to be about  
// three times slower than those involving only multiplication and division.
//
{
    n_cycles = 0 ;
    max_diff = 2.0 * tol ; // Dummy value: replaced before use
	while( (n_cycles == 0) || ((max_diff > tol) && (n_cycles < max_cycles)) ) {
        // Eqn. (4.43)
        Sfg = C ;
        Sfg.each_row() -= g.t() ;
        Sfg.each_col() -= f ;
        f += soft_mins_across_rows(col_tmp, epsilon, Sfg) + eps_log_a ;
        
        // Eqn. (4.44)
        Sfg = C ;
        Sfg.each_row() -= g.t() ;
        Sfg.each_col() -= f ;
        g += soft_mins_down_cols(row_tmp, epsilon, Sfg) + eps_log_b ;
               
        // Recompute the solution and the convergence score
        // every cycles.per.check cycles.
        n_cycles++ ;
        if( n_cycles % cycles_per_check == 0 ) {
			log_P = neg_C_over_eps ;
			log_P.each_row() += g.t() / epsilon ;
			log_P.each_col() += f / epsilon ;
			P = arma::exp(log_P) ;
			arma::colvec row_margin_abs_diff = arma::abs(arma::sum(P, 1) - a) ;
			max_diff = row_margin_abs_diff.max() ;
			converged = max_diff < tol ;
			if( printProgress ) {
				Rcout << n_cycles << ": max_diff=" << max_diff ;
				Rcout << ", sum(P)=" << arma::accu(P) << std::endl ;
			}
		}
    }
}