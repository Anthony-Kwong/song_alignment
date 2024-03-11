#' og_order function
#' 
#' Current version of the align function changes order of the sequences. This function will look at
#' the alignment matrix A and the original sequences S; for every song in S, it will find the corresponding
#' row in A. Return a new matrix A_star, consisting of the rows in A but in the order in S.
#'
#' @param align_mat : character matrix
#' @param song_seqs : string vector
#'
#' @return Character matrix A_star with rows of align_mat, in the order of the sequences in song_seqs.
#' @export
#'
#' @examples
og_order <- function(align_mat, song_seqs) {
  #turn rows of align_mat back into strings, removing the gaps
  string_vectors <- apply(align_mat, 1, function(r){
    x = paste(r, collapse = "")
    cleaned_x = gsub("-", "", r)
    paste(cleaned_x, collapse = "")
  })
  
  #match string_vec with the songs in song_seqs
  indices = sapply(song_seqs, function(r){
    #in the case where 2 sequences are the same, we only need to retrieve the first index
    which(string_vectors == r)[1]
  })
  indices = as.vector(indices)
  A_star = align_mat[indices,]
  
  return(A_star)
}

#tests 

song_seqs = c("ACCBC", "ABBA", "ABBA")
A = t(matrix(c("A","B","B","A","-", "A", "C", "C", "B", "C", "A","B","B","A","-"), ncol = 3))
output = og_order(align_mat = A, song_seqs = song_seqs)
testthat::expect_equal(A[c(2,1,1),], output)
