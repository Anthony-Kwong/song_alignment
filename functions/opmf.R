#' opmf function
#' 
#' Computes the pmf of the offset values in the Gibbs aligner.
#'
#' @param o : numeric integer for the offset
#' @param P : Propensity matrix
#' @param t_star : character string for the motif. Each element is a letter.
#' @param shift : n_star-w
#'
#' @return Numeric. Probability mass at offset o. 
#' @export
#'
#' @examples
opmf <- function(o, P, t_star, shift){
  #compute numerator----
  top = sapply(seq(w), function(j){
    l = which(rownames(P) == t_star[o+j])
    P[l,j]
  })
  top_prod = prod(top)
  #compute denominator----
  bot = lapply(1:(shift+1), function(i) i)
  for(i in 0:shift){
    bot[[i+1]] = sapply(seq(w), function(j){
      l = which(rownames(P) == t_star[i+j])
      P[l,j]
    })
  }
  prod_bot = lapply(bot, function(x){prod(x)})
  sum_bot = sum(as.numeric(prod_bot))
  #compute final ratio
  res = top_prod/sum_bot
  return(res)
}

#write test....

#construct a biggish tstar without too much shift

#check probs sum to 1