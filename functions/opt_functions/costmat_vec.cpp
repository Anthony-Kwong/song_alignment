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

//takes 2 numeric vectors and generates a cost matrix for optimal transport problem

//input: Numeric vectors x,y denoting the ticks of the bins
//output: Numeric cost matrix of Wassertein distances

// [[Rcpp::export]]
NumericMatrix costmat_vec(NumericVector x, NumericVector y) {
  int xlen = x.length();
  int ylen = y.length();
  
  //intialise cost matrix
  NumericMatrix C(xlen, ylen);
  //loop rows
  for(int i = 0; i<xlen; i++){
    //loop cols
    for(int j = 0; j<ylen; j++){
      C(i, j) = std::pow(x[i] - y[j], 2);
    }
  }
  return C;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#test
x = seq(4)
output = costmat_vec(x,x)

dims = length(x)
ans = matrix(ncol=dims, nrow=dims)
for(i in 1:dims){
  for(j in 1:dims){
    ans[i,j] = (x[i]-x[j])^2
  }
}

testthat::expect_equal(ans,output)
*/
