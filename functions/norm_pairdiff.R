#' norm_pairdiff function
#' 
#' Takes a set of aligned songs from 2 birds and computes the normalised pairwise differences between the 2 birds'
#' songs. The normalisation term is the number of pairwise comparisons made. 
#'
#'
#' @param s1 : Character matrix showing the aligned songs of bird1.
#' @param s2 : Character matrix showing the aligned songs of bird2.
#'
#' @return Numeric. Number of pairwise differences between the 2 birds; normalized. 
#' @export
#'
#' @examples

source("./functions/string_differences.R")
norm_pairdiff <- function(s1,s2){
  #add dummy index variables to get all combinations of pairs for comparison----
  n1 = seq(nrow(s1))
  n2 = seq(nrow(s2))
  
  #get all unique pairwise combinations----
  pair_combs = expand.grid(n1, n2)
  pairs = nrow(pair_combs)
  #perform comparison for every pair----
  
  #preallocate vector of differences in individual pairs
  diff = rep(NA, pairs)
  for(i in 1:nrow(pair_combs)){
    pair = pair_combs[i,]
    diff[i] = string_differences(s1[pair$Var1,], s2[pair$Var2,])
  }
  
  #return normalised result
  ans = mean(diff)
  return(ans)
}

#tests
a = rbind(c("A","B","B","A"), c("C","B","B","A"))
b = rbind(c("A","D","-","A"), c("A","D","B","A"))

total_diff = 8
norm_diff = 8/4
output = norm_pairdiff(a,b)
testthat::expect_equal(output, norm_diff)