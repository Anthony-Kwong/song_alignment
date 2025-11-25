#' string_slice function
#' 
#' Takes an input string and returns a substring using start and end indices. 
#'
#' @param s : A character string
#' @param start : Index of starting character of substring. 
#' @param end : Index of ending character for substring.
#'
#' @return : substring of s, bounded by start and end. 
#' @export
#'
#' @examples s = "ABBACC", string_slice(s, 2, 4)
string_slice <- function(s, start, end){
  #check input
  if(is.character(s)==FALSE){
    stop("input string s must be character.")
  }
  if(start > end){
    stop("start index must be before end index")
  }
  
  #split string to character vector
  S = unlist(strsplit(s, split=""))
  sub_S = substring = S[start:end]
  return(sub_S)
}

testthat::test_that("string_slice works",{
  s = "ABBACC"
  sub_S = string_slice(s, 1, 4)
  ans = c("A", "B", "B", "A")
  testthat::expect_equal(sub_S, ans)
})

testthat::test_that("string_slice works",{
  s = "CGABDKASDFCX"
  sub_S = string_slice(s, 4, 7)
  ans = c("B", "D", "K", "A")
  testthat::expect_equal(sub_S, ans)
})

