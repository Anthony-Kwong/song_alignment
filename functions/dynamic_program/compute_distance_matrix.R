#' Compute distance matrix function
#' 
#' Compute the distance matrix of a sequence set, using a 
#'
#' @param S : A list of character vectors. 
#' @param gap : gap penalty score
#' @param mismatch : mismatch penalty score
#' @param match : match score
#' @param method: dynamic programming method. Either needleman-wunsch(NW) or smith-waterman(SW). Default NW. 
#'
#' @return Numeric distance matrix, where the elements are the pairwise alignment scores
#' @export
#'
#' @examples
compute_distance_matrix <- function(S, gap, mismatch, match, method = "NW") {
  #check input method
  valid_methods = c("NW","SW")
  if(!(method %in% valid_methods)){
    stop("invalid method. method must be either NW or SW")
  }
  #get nseqs
  n <- length(S)
  #call distance mat
  D <- matrix(-Inf, n, n)
  
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      if(method == "SW"){
        #smith-waterman method
        aln = text.alignment::smith_waterman(a = S[i], b = S[j], match = match, mismatch = mismatch, gap = gap)
        score = aln$sw
      } else {
        #needleman-wunsch method (default)
        aln = needleman(seq1 = S[i], seq2 = S[j], gap = gap, mismatch = mismatch, match = match)
        score = aln$score
      }
      D[i,j] <- score
    }
  }
  return(D)
}

testthat::test_that("",{
  S = c("ABBACC", "CCABBA", "CABBA")
  output = compute_distance_matrix(S, match = 2, gap = -1, mismatch = -1,method = "NW")
  
  #checked scores manually
  D <- matrix(-Inf, n, n)
  D[1,2] = 4
  D[1,3] = 5
  D[2,3] = 9
  
  testthat::expect_equal(output, D)
})


