#' Progressive alignment function
#' 
#' Progressive alignment using pairwise dynamic programming to make multiple sequence alignments.
#'
#' @param S : A list of character vectors. 
#' @param match : match score. Default 2. 
#' @param mismatch : mismatch score. Default -1. 
#' @param gap : gap penalty. Default -1. 
#'
#' @return Character matrix for multiple sequence alignment of sequences S. 
#' @export
#'
#' @examples S = c("ABBACC", "CCABBA", "CABBA")
#' progressive_align(S)

source("./functions/needleman.R")
source("./functions/dynamic_program/compute_distance_matrix.R")
progressive_align <- function(S, match = 2, mismatch = -1, gap = -1){
  
  
  # Initialisation. No aln yet. Store aligned sequences in aln_mat.
  aln_mat <- NULL
  sequences <- S
  
  #counter for tracking
  counter = 0
  
  #iterate until we have added all sequences into the alignment
  while(length(sequences) > 0) {
    
    if(is.null(aln_mat)) {
      # First iteration — find best pair
      D <- compute_distance_matrix(sequences, gap, mismatch, match, method = method)
      #get index of the highest scoring pair, i j
      ij <- which(D == max(D), arr.ind=TRUE)[1,]
      #get indices of seqs i,j
      i <- ij[1]; j <- ij[2]
      #pairwise align i,j
      
      #needleman-wunsch default
      aln <- needleman(sequences[i], sequences[j],
                       gap = gap, mismatch = mismatch, match = match)
      aln_mat <- aln$aligned_seqs
      
      # remove used sequences
      keep <- setdiff(1:length(sequences), c(i,j))
      sequences <- sequences[keep]
      
    } else if(length(sequences) > 0) {
      # Align next sequence to current consensus
      
      # Current consensus string: ---- 
      cons <- apply(aln_mat, 2, function(col) {
        # remove gaps
        col <- col[col != "-"]
        # if column completely gaps, this should not trigger
        if (length(col) == 0) return("-")
        # count frequencies
        tab <- table(col)
        maxfreq <- max(tab)
        # letters with max frequency, (allowing for a tie)
        winners <- names(tab)[tab == maxfreq]
        # if tie → choose one randomly
        sample(winners, 1)
      }) |> paste(collapse="")
      
      # Find the next highest scoring sequence to pair with the consensus ----
      
      #initialization
      best <- NULL
      bestScore <- -Inf
      #index for best sequence
      best_i <- NULL
      
      #loop along remaining sequences
      for(i in seq_along(sequences)) {
        
        #default needleman method
        aln <- needleman(seq1 = cons, seq2 = sequences[i], gap = gap, mismatch = mismatch, match = match)
        score = aln$score
        
        if(score > bestScore) {
          bestScore <- score
          best_i <- i
          best <- aln
        }
      }
      
      # Now align this sequence to the existing MSA
      
      # Extract aligned consensus and sequence
      cons_aln <- best$aligned_seqs[1,]
      seq_aln  <- best$aligned_seqs[2,]
      
      # We must insert gaps into aln_mat according to cons_aln, incase the pairwise alignment adde a gap 
      gap_positions <- which(cons_aln == "-")
      
      if(length(gap_positions) > 0) {
        #initialise new matrix
        new_mat <- matrix("-", nrow=nrow(aln_mat), ncol=length(cons_aln))
        idx <- setdiff(seq_along(cons_aln), gap_positions)
        #This shifts all original columns into the appropriate positions, leaving gaps in the correct new gap locations.
        new_mat[, idx] <- aln_mat
        #update aln_mat
        aln_mat <- new_mat
      }
      
      # Append the new aligned sequence
      aln_mat <- rbind(aln_mat, seq_aln)
      print(paste("added a seq",counter))
      counter = counter + 1
      
      # Remove it from the pool
      sequences <- sequences[-best_i]
    }
  }
  return(aln_mat)
}

testthat::test_that("",{
  #we will only check that the output is valid
  
  input_sequences = c("ABBACC", "CCABBAD", "CABBCA")
  output = progressive_align(S)
  
  # Remove gaps from every row in the aligned matrix
  stripped <- apply(output, 1, function(row) paste(row[row != "-"], collapse=""))
  
  # Check that all stripped sequences are in original set
  all_in_original <- all(stripped %in% input_sequences)
  testthat::expect_equal(all_in_original, T)
  
  # Check that counts match — no repeats, no losses
  same_set <- setequal(stripped, input_sequences)
  testthat::expect_equal(same_set, T)
})

# #testing output is valid
# test_alignment <- function(aligned_matrix, input_sequences) {
#   
#   # Remove gaps from every row in the aligned matrix
#   stripped <- apply(aligned_matrix, 1, function(row) paste(row[row != "-"], collapse=""))
#   
#   # Check that all stripped sequences are in original set
#   all_in_original <- all(stripped %in% input_sequences)
#   
#   # Check that counts match — no repeats, no losses
#   same_set <- setequal(stripped, input_sequences)
#   
#   # Output results
#   list(
#     stripped_sequences = stripped,
#     are_valid_sequences = all_in_original,
#     set_matches_original = same_set
#   )
# }

# a = bird_songs %>%
#   dplyr::filter(Line == "Pink")
# asongs = a$note.seq
# 
# k = progressive_align(asongs[1:10])
