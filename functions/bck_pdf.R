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
  #get alphabet
  alphabet = names(prop_tab)
  prop_mat = matrix(prop_tab)
  rownames(prop_mat) = alphabet
  return(prop_mat)
}

testthat::test_that("",{
  data = rbind(c("A","B","B","A"), c("C","C","C","A"),c("A","B","C","C"))
  output = bck_pdf(data)
  
  ans = matrix(c(4/12,3/12,5/12))
  rownames(ans) = c("A","B","C")
  testthat::expect_equal(output, ans)
})


