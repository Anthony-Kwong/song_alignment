#' Progressive alignment function (deprecated)
#' 
#' This version was attemmpting to get SW to work as well. 
#' 
#' Progressive alignment using pairwise dynamic programming to make multiple sequence alignments.
#'
#' @param S : A list of character vectors. 
#' @param match : match score. Default 2. 
#' @param mismatch : mismatch score. Default -1. 
#' @param gap : gap penalty. Default -1. 
#' @param method : Method of dynamic programming. Options are Needleman-Wunsch(NW) and Smith-Waterman(SW). Default NW. 
#' SW nolonger suppored. 
#'
#' @return Character matrix for multiple sequence alignment of sequences S. 
#' @export
#'
#' @examples S = c("ABBACC", "CCABBA", "CABBA")
#' progressive_align(S)

source("./functions/needleman.R")
source("./functions/dynamic_program/compute_distance_matrix.R")
progressive_align <- function(S, match = 2, mismatch = -1, gap = -1, method = "NW"){
  
  
  # Initialisation. No aln yet. Store aligned sequences in aln_mat.
  aln_mat <- NULL
  sequences <- S
  
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
      
      #smith-waterman method
      if(method == "SW"){
        aln <- text.alignment::smith_waterman(a = sequences[i], b = sequences[j],
                                              gap = gap, mismatch = mismatch, match = match)
        #retrieve aligned seqs
        a = aln$a$tokens
        b = aln$b$tokens
        aln_mat = rbind(a,b)
      } else {
        #needleman-wunsch default
        aln <- needleman(sequences[i], sequences[j],
                         gap = gap, mismatch = mismatch, match = match)
        aln_mat <- aln$aligned_seqs
      }
      
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
        if(method == "SW"){
          #smith-waterman
          aln = text.alignment::smith_waterman(a = cons, b = sequences[i],
                                               gap = gap, mismatch = mismatch, match = match)
          score = aln$sw
          aln = rbind(aln$a$tokens, aln$b$tokens)
        } else {
          #default needleman method
          aln <- needleman(cons, sequences[i], gap = gap, mismatch = mismatch, match = match)
          score = aln$score
        }
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
      
      # Remove it from the pool
      sequences <- sequences[-best_i]
    }
  }
  return(aln_mat)
}