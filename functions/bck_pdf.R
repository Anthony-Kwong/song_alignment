#' bck_pdf function
#' 
#' Takes a set of characters strings and returns the proportion of every occurring letter.
#'
#' @param data: A list of character strings. 
#'
#' @return: A table containing the proportions of every occurring letter. 
#' @export
#'
#' @examples data = c("CCABBACCCC", "DDABBADDDD", "EEABBAEEEE")
#' bck_pdf(data)
bck_pdf <- function(data){
  #combine all letters into one big string
  seqs = unlist(strsplit(data, split = ""))
  letter_counts = table(seqs)
  prop_tab = prop.table(letter_counts)
  return(prop_tab)
}

library(testthat)
data = c("CCABBACCCC", "DDABBADDDD", "EEABBAEEEE")
b = bck_pdf(data)
expect_equal(as.numeric(b), rep(0.2,5))

