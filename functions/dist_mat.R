#' dist_mat function
#' 
#' Takes a table of bird pairs and their distances; generates a distance matrix using those values. 
#'
#' @param tab : A 3 column dataframe with columns (bird1, bird2, distance)
#'
#' @return : A Numeric distance matrix for the birds.
#' @export
#'
#' @examples See tests
dist_mat <- function(tab){
  #get the names of all the birds
  birds = unique(c(tab$V1,tab$V2))
  #initialize matrix
  matrix <- matrix(0, nrow = length(birds), ncol = length(birds),
                   dimnames = list(birds, birds))
  k = colnames(matrix)
  #add entries to matrix(we will do upper diagonal)
  for(i in 1:nrow(tab)){
    sam = tab[i,]
    i = which(k==sam$V1)
    j = which(k==sam$V2)
    matrix[i,j] = sam$D
  }
  
  #matrix is symmetric
  ans = matrix + t(matrix)
  return(ans)
}

#tests