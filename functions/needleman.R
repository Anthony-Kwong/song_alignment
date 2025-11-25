# References: https://de.wikipedia.org/wiki/Needleman-Wunsch-Algorithmus

#Function written by Julius Kittler

#Sourced from https://gist.github.com/juliuskittler/ed53696ac1e590b413aac2dddf0457f6. 

needleman = function(seq1, seq2, gap, mismatch, match){
  
  # Stop conditions
  stopifnot(gap <= 0) # check if penalty negative
  stopifnot(mismatch <= 0)  # check if penalty negative
  stopifnot(match >= 0)  # check if score positive
  
  # Initialize col and rownames for matrices
  len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
  seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
  seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
  
  # Initialize matrix M (for scores)
  M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(M) = c("-", seq1) # assign seq chars to matrix names
  colnames(M) = c("-", seq2) # assign seq chars to matrix names
  M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
  M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
  
  # Initialize matrix D (for directions)
  D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(D) = c("-", seq1) # assign seq chars to matrix names
  colnames(D) = c("-", seq2) # assign seq chars to matrix names
  D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
  D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
  type = c("dia", "hor", "ver") # Lookup vector
  
  # Compute scores and save moves
  for (i in 2:(len1 + 1)){# for every (initially zero) row
    for (j in 2:(len2 + 1)){# for every (initially zero) col
      hor = M[i, j - 1] + gap # horizontal move = gap for seq1
      ver = M[i - 1, j] + gap # vertical move = gap for seq2 
      dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch) 
                   M[i - 1, j - 1] + match, 
                   M[i - 1, j - 1] + mismatch)
      M[i, j] = max(dia, hor, ver) # Save current (best) score in M
      D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
    }
  } 
  
  # Backtracing
  align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
  
  #reset backtrace starting point
  # i = len1 + 1
  # j = len2 + 1
  
  while(i > 1 || j > 1){
    
    if(D[i, j] == "dia") {
      align1 = c(rownames(M)[i], align1)
      align2 = c(colnames(M)[j], align2)
      j = j - 1; i = i - 1  # update indices
    } else if (D[i, j] == "ver") {
      align1 = c(rownames(M)[i], align1)
      align2 = c("-", align2) # vertical movement = gap for seq2
      i = i - 1 # update indices
    } else if (D[i, j] == "hor") {
      align1 = c("-", align1) # horizontal movement = gap for seq1
      align2 = c(colnames(M)[j], align2) 
      j = j - 1 # update indices
    } 
    
  }
  
  # Prepare output
  return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
              score = M[nrow(M), ncol(M)], score_matrix = M, movement_matrix = D))
  
}

# Test case 1: Wiki example (https://de.wikipedia.org/wiki/Needleman-Wunsch-Algorithmus)
testthat::test_that("needlemanreturns correct alignment", {
  solution = needleman("ACGTC", "AGTC", gap = -1, mismatch = -1, match = 0)$aligned_seqs
  expected = matrix(c("A", "C", "G", "T", "C", 
                      "A", "-", "G", "T", "C"), byrow = TRUE, nrow =2)
  testthat::expect_equal(solution, expected)
})

# Test case 2: Problem 2.9 (M. Borodovsky, S. Ekisheva., 2006, 
# Problems and Solutions in Biological Sequence Anal- ysis, Cambridge University Press)
testthat::test_that("needlemanreturns correct alignment", {
  solution = needleman("GAATTC", "GATTA", gap = -2, mismatch = -1, match = 2)$aligned_seqs
  expected = matrix(c("G", "A", "A", "T", "T", "C",
                      "G", "-", "A", "T", "T", "A"), byrow = TRUE, nrow =2)
  testthat::expect_equal(solution, expected)
})

# Test case 3: Fig. 5-7 (J. Momand, A. McCurdy., 2017, Concepts in Bioinformatics  
# and Genomics, Oxford Uni- versity Press.)
testthat::test_that("needlemanreturns correct alignment", {
  solution = needleman("ADCDNRCKCRWP", "AWCNDRQCLCRP", 
                       gap = 0, mismatch = 0, match = 1)$aligned_seqs
  expected = matrix(c("A","D","C","D","N","-","R","-","C","K","C","R","W","P",
                      "A","W","C","-","N","D","R","Q","C","L","C","R","-","P"), 
                    byrow = TRUE, nrow =2)
  testthat::expect_equal(solution, expected)
})


# Test case 4: Ch. 5, Problem 5 MODIFIED (!) (J. Momand, A. McCurdy., 2017, Concepts in 
# Bioinformatics and Genomics, Oxford Uni- versity Press.)
testthat::test_that("needlemanreturns correct alignment", {
  solution = needleman("ATAGC", "ATATGA", gap = 0, mismatch = 0, match = 1)$aligned_seqs
  expected = matrix(c("A","T","A","-","G","C" ,
                      "A","T","A","T","G","A"), byrow = TRUE, nrow =2)
  testthat::expect_equal(solution, expected)
})

