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
NumericMatrix costmat_C(NumericVector pf, NumericVector pt, NumericVector qf, NumericVector qt) {
  
  int plen = pf.length();
  int qlen = qf.length();
  
  std::cout << plen << std::endl; 
  std::cout << qlen << std::endl; 
  
  
  //intialise cost matrix
  NumericMatrix C(qlen, plen);
  //loop rows
  for(int i = 0; i<qlen; i++){
    //loop rows
    for(int j = 0; j<plen; j++){
      C(i, j) = std::pow(qf[i] - pf[j], 2) + std::pow(qt[i] - pt[j], 2);
    }
  }
    //(Q$freq.ticks - P$freq.ticks)^k + a^(k)*(Q$time.ticks - P$time.ticks)^k
  return C;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#D = costmat_C(p_vals$freq.ticks, p_vals$time.ticks, q_vals$freq.ticks, q_vals$time.ticks)

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
