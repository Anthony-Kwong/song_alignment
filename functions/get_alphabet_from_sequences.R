#' get_alphabet_from_sequences functions
#' 
#' Take a character vector and returns all unique letters.
#'
#' @param sequences : character vector
#'
#' @return : character vector of all unique letters in sequences
#' @export
#'
#' @examples get_alphabet_from_sequences(c("ABBACC", "CCABBA", "CABBA"))
get_alphabet_from_sequences <- function(sequences) {
  # Split each sequence into characters
  split_letters <- strsplit(sequences, split = "")
  
  # Combine all letters into one vector
  all_letters <- unlist(split_letters)
  
  # Remove any empty strings just in case
  all_letters <- all_letters[all_letters != ""]
  
  # Get sorted unique set of letters
  alphabet <- sort(unique(all_letters))
  
  return(alphabet)
}

testthat::test_that("",{
  output = get_alphabet_from_sequences(c("ABBACC", "CCABBA", "CABBA"))
  ans = c("A","B","C")
  testthat::expect_equal(output, ans)
})
