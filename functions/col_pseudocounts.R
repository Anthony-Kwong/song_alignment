#' col_pseudocounts function
#' 
#' Takes a vector of letters and returns the distribution using their frequencies via the pseudocounts method.
#'
#' @param let_vec : Character vector of letters. (Columns from the alignment matrix)
#' @param letters : Character vector. Set of letters which appear in the overall alignment.
#' @param pcount : Numeric. Pseudocount. 
#' @param pro: Logical (default TRUE). Set to FALSE to return counts rather than proportions.
#'
#' @return Probability distribution for the input vector.
#' @export
#'
#' @examples See tests.
col_pseudocounts <- function(let_vec, letters, pcount, pro = TRUE){
  #remove the gaps
  let_vec <- let_vec[let_vec != "-"]
  
  #calculate frequencies in let_vec
  counts <- table(let_vec)
  #add pseudocounts
  for (letter in letters) {
    if (!(letter %in% names(counts))) {
      counts[letter] <- 0  # Add zero count if the letter is not in V
    }
    counts[letter] <- counts[letter] + pcount
  }
  
  sorted_counts <- counts[order(names(counts))]
  
  #return counts case
  if(pro){
    # Calculate total count including pseudocounts
    total_count <- sum(sorted_counts)
    # Compute proportions
    proportions <- sorted_counts / total_count
    return(proportions)
  } else {
    return(as.numeric(sorted_counts))
  }
}

#testing

# Sample vector V
V <- c("A", "D", "A", "D", "A", "A", "A", "D", "D")
# Letters to consider
l <- LETTERS[1:4]
output = col_pseudocounts(let_vec = V, letters = l, pcount = 1)
p_freq = c(5,0,0,4) + 1
ans = p_freq/sum(p_freq)
testthat::expect_equal(as.numeric(output), ans)

#test 2
V = c("B", "B", "A", "B", "B", "A", "A", "D", "D")
l <- LETTERS[1:5]
output = col_pseudocounts(let_vec = V, letters = l, pcount = 1)
p_freq = c(3,4,0,2,0) + 1
ans = p_freq/sum(p_freq)
testthat::expect_equal(as.numeric(output), ans)

#test 3

V = c("B", "A", "A", "B", "C", "A", "A", "D", "D")
l <- LETTERS[1:5]
output = col_pseudocounts(let_vec = V, letters = l, pcount = 1, pro = FALSE)
ans = c(4,2,1,2,0) + 1
testthat::expect_equal(output, ans)
