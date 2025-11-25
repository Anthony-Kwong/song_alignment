#' string_split function
#' 
#' Takes in a vector of character strings and returns a list, with every string split into separate characters.
#'
#' @param seqs : Vector of character strings
#'
#' @return : List of character vectors, with every string split into its constituent characters.
#' @export
#'
#' @examples S = c("ABCAC", "ABCA")
#' string_split(S)
string_split <- function(seqs){
  #get sequences in split strings
  seqs_split = lapply(seqs, function(s){strsplit(s, "")[[1]]})
  return(seqs_split)
}

testthat::test_that("",{
  S = c("ABCAC", "ABCA")
  output = string_split(S)
  
  ans = list(c("A","B","C","A","C"), c("A","B","C","A"))
  testthat::expect_equal(output,ans)
})
