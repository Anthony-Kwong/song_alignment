#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//input: 
// pf : NumericVector for the frequency ticks in source distribution
// pt : NumericVector for the time ticks in source distribution
// qf : NumericVector for the frequency ticks in destination distribution
// qt : NumericVector for the time ticks in destination distribution

// [[Rcpp::export]]
NumericMatrix costmat_C(NumericVector pf, NumericVector pt, NumericVector qf, NumericVector qt, double a=1.0, double k=2.0 ) ;
NumericMatrix costmat_C(NumericVector pf, NumericVector pt, NumericVector qf, NumericVector qt, double a, double k ) {
  
  int plen = pf.length();
  int qlen = qf.length();
  
  Rcout << plen << ", " << qlen ;
  Rcout << ": "<< (plen/256)*(qlen/256) << std::endl ; 
  
  try{
    //intialise cost matrix
    NumericMatrix C(qlen, plen);
    //loop rows
    for(int i = 0; i<qlen; i++){
      //loop cols
      for(int j = 0; j<plen; j++){
        C(i, j) = std::pow(qf[i] - pf[j], k) + std::pow(a*(qt[i] - pt[j]), k);
      }
    }
    //(Q$freq.ticks - P$freq.ticks)^k + a^(k)*(Q$time.ticks - P$time.ticks)^k
    return C;
  }
  catch(std::exception &ex) {
    return(R_NilValue) ;
  }
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#run tests

pf = seq(1,22, by = 5)
pt = seq(5)
qf = seq(1,22, by = 5)
qt = seq(5)

output = costmat_C(pf, pt, qf, qt)

#loop cols
for(i in 1:5){
  ans = (pf[i]-qf)^2 + (pt[i]-qt)^2
  testthat::expect_equal(ans, output[i,])
}

*/
