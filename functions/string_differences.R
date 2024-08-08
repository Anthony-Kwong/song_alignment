#' calculate_differences function
#' 
#' Take 2 character vectors and return the number of pairwise differences. Elements can include letters and "-" (gaps).
#'
#' @param vector1 : Character vector.
#' @param vector2 : Character vector. 
#'
#' @return Numeric. Number of pairwise differences between the 2 character strings. 
#' @export
#'
#' @examples see tests
string_differences <- function(vector1, vector2) {
  # Check if vectors have the same length
  if (length(vector1) != length(vector2)) {
    stop("Vectors must have the same length")
  }
  
  # Initialize counter for differences
  differences <- 0
  
  # Iterate over corresponding elements of vectors and count differences
  for (i in 1:length(vector1)) {
    if (vector1[i] != vector2[i]) {
      differences <- differences + 1
    }
  }
  
  # Return the number of differences
  return(differences)
}

#tests----
v1 <- c("A", "N", "D", "C")
v2 <- c("A", "N", "B", "C")
output = string_differences(v1,v2)
testthat::expect_equal(output, 1)

#check distance of 0
output = string_differences(v1,v1)
testthat::expect_equal(output, 0)

v1 <- c("A", "-", "B", "-", "D")
v2 <- c("A", "-", "B", "C", "D")
output = string_differences(v1,v2)
testthat::expect_equal(output, 1)

