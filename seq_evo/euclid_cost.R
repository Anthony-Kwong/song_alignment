#' compute euclidean cost matrix for 2 vectors
#'
#' @param x : Numeric vector
#' @param y : Numeric vector
#'
#' @return : Cost matrix C where Cij is the euclid distance between Xi and Yj
#' @export
#'
#' @examples
euclid_cost <- function(x, y){
  C = matrix(nrow = length(x), ncol = length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      C[i,j] = x[i]-y[j]
    }
  }
  return(C)
}
