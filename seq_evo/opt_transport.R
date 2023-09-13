#infer sequence evolution using optimal transport

#' Discretised optimal transport function
#' 
#' Compute optimal transportation plan between 2 pmfs.
#'
#' @param pmf1 : Probability mass function expressed as a vector. 
#' @param pmf2 : Probability mass function expressed as a vector. 
#' @param cost : Cost function expressed as a numeric matrix. 
#'
#' @return Optimal transport plan as a numeric matrix
#' @export
#'
#' @examples
opt_transport <- function(pmf1, pmf2, cost){
  obj.fn = as.vector(t(cost))
  
  # Assemble the matrix of constraint coefficients. We use a sparse representation
  #Construct constraints matrix----
  n = length(pmf1)
  #constraints T1 = p
  #initialize for T1 = p
  const.mat1 = matrix(0, ncol = n*n, nrow = n)
  #loop across rows
  for(i in 1:(nrow(const.mat1))){
    print(i)
    #const.mat1[i,][(i+n*(i-1)):(i+n*i-1)] = 1
    const.mat1[i,][(1+n*(i-1)):(n*i)] = 1
  }

  #constraints T1 = q
  const.mat2 = matrix(0, ncol = n*n, nrow = n)
  for(j in 1:(nrow(const.mat1))){
    index = seq(from=j, to = n*n, by=n)
    const.mat2[j,][index] = 1
  }
  #combine the constraint matrices
  const.mat = rbind(const.mat1,const.mat2)
  
  #set constraint directions
  const.dir <- rep("<=",2*n)
  
  #set rhs of constraints
  const.rhs = c(pmf1, pmf2)
  
  lp.solution <- lpSolve::lp("min", obj.fn, const.mat, 
                    const.dir, const.rhs, compute.sens=TRUE)
  
  return(lp.solution)
}

# # consisting of (row, col, value) triples, that lpSolve() accepts.
# colsum.constraint.sparse.mat <- cbind( 
#   rep(1:n, each=n), # constraint number
#   1:(n*n), # index of variable
#   rep(1, n*n)
# )