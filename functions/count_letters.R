#' count_letters
#' 
#' Count every letter in a character vector, given some set of letters.
#'
#' @param letters : character vector
#' @param alphabet : set of letters
#'
#' @return Numeric vector containing the count of every letter in alphabet.
#' @export
#'
#' @examples letters <- c("A", "B", "A", "C", "A")
#' alphabet <- c("A", "B", "C", "D")
#' count_letters(letters,alphabet)

count_letters <- function(letters, alphabet) {
  # Ensure input is character
  letters <- as.character(letters)
  alphabet <- as.character(alphabet)
  
  # Create a frequency table of observed letters
  freq <- table(factor(letters, levels = alphabet))
  
  # Convert to numeric named vector for easy handling
  counts <- as.numeric(freq)
  names(counts) <- alphabet
  
  return(counts)
}

testthat::test_that("",{
  letters <- c("A", "B", "A", "C", "A")
  alphabet <- c("A", "B", "C", "D")
  output = count_letters(letters,alphabet)
  
  ans = c(3,1,1,0)
  testthat::expect_equal(ans, as.numeric(output))
})


testthat::test_that("",{
  letters <- c("A", "B", "C", "C", "A", "E", "E")
  alphabet <- c("A", "B", "C", "D", "E")
  output = count_letters(letters,alphabet)
  
  ans = c(2,1,2,0,2)
  testthat::expect_equal(ans, as.numeric(output))
})


