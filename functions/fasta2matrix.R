#' fasta2matrix function
#' 
#' Turns a fasta file from a directory and ouputs alignment matrix
#'
#' @param file_path : file path to the fasta file
#'
#' @return : Character matrix 
#' @export
#'
#' @examples see tests
fasta2matrix <- function(filepath) {
  # Read the entire file
  fasta_lines <- readLines(filepath)
  
  # Initialise containers
  seq_names <- character()
  sequences <- list()
  
  # Parse the FASTA file
  current_name <- NULL
  current_sequence <- ""
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      # If there's an ongoing sequence, save it
      if (!is.null(current_name)) {
        seq_names <- c(seq_names, current_name)
        sequences <- c(sequences, list(current_sequence))
      }
      # Start a new sequence
      current_name <- sub(">", "", line)
      current_sequence <- ""
    } else {
      # Concatenate sequence lines
      current_sequence <- paste0(current_sequence, line)
    }
  }
  # Save the last sequence
  if (!is.null(current_name)) {
    seq_names <- c(seq_names, current_name)
    sequences <- c(sequences, list(current_sequence))
  }
  
  # Convert sequences to a matrix
  max_length <- max(nchar(sequences))
  sequence_matrix <- do.call(rbind, lapply(sequences, function(seq) {
    strsplit(seq, "")[[1]]
  }))
  rownames(sequence_matrix) <- seq_names
  
  return(sequence_matrix)
}

# Small example for testing (by eye)
test_fasta <- tempfile(fileext = ".txt")
writeLines(c(
  ">seq1", "AA--AAAAA--FFGAFAFFFAAAAAAAFGAFAFF-FA-AAAAA-FG",
  ">seq2", "AA--AAAAF--FFGAFAFFFAAAAAAAFGAFAFF-FA-AAAAA-FG",
  ">seq3", "A---AAAAF--FFGAFAFFFAA-AAAAFGAFAFF-FA-AAAAA-FG"
), test_fasta)

# Test the function
fasta_matrix <- fasta2matrix(test_fasta)
print(fasta_matrix)

#proper computational test
# Create a test FASTA file
test_fasta <- tempfile(fileext = ".txt")
writeLines(c(
  ">seq1", "AA--A",
  ">seq2", "AA--G",
  ">seq3", "AA--T"
), test_fasta)

# Run the function
result <- fasta2matrix(test_fasta)

# Expected matrix
expected <- matrix(
  c("A", "A", "-", "-", "A",
    "A", "A", "-", "-", "G",
    "A", "A", "-", "-", "T"),
  nrow = 3, byrow = TRUE
)
rownames(expected) <- c("seq1", "seq2", "seq3")

# Test if the result matches the expected output
testthat::expect_equal(result, expected)



