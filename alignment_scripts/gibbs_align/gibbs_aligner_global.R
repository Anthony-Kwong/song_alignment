#' Gibbs aligner for arbitrary string sequences
#'
#' This function performs a simple Gibbs-sampling motif alignment on a vector
#' of strings. It samples motif start positions of fixed width `w`, then
#' returns a full alignment that preserves the unaligned prefixes and suffixes
#' by padding with `-` so the sampled motif columns line up.
#'
#' @param S Vector of strings to align.
#' @param w Numeric integer. Width of alignment / motif.
#' @param iter Numeric scalar. Number of Gibbs sweeps to carry out.
#'
#' @return A character matrix of the final alignment.
#'
#' @details
#' The output is a rectangular character matrix in which all rows have the same
#' number of columns. The aligned width-`w` windows are placed in the same
#' columns, while any unaligned characters at the left or right ends are kept
#' and shifted accordingly. Empty positions are filled with `-`.
#'
#' For example, if the chosen aligned window is `ABBA` for
#' `c("XABBA", "XXABBA", "ABBAX")`, one valid final alignment is:
#'
#' `-XABBA-`
#' `XXABBA-`
#' `--ABBAX`
#'
#' @examples
#' set.seed(1)
#' x <- gibbs_aligner_full(c("XABBA", "XXABBA", "ABBAX"), w = 4, iter = 200)
#' apply(x, 1, paste0, collapse = "")
#'
#' @export
gibbs_aligner_global <- function(S, w, iter) {
  # Validate the sequence input.
  if (!is.character(S) || length(S) < 2L) {
    stop("S must be a character vector of length >= 2.")
  }
  
  # Validate motif/alignment width.
  if (length(w) != 1L || !is.finite(w) || w < 1 || w != as.integer(w)) {
    stop("w must be a positive integer scalar.")
  }
  
  # Validate number of Gibbs sweeps.
  if (length(iter) != 1L || !is.finite(iter) || iter < 1 || iter != as.integer(iter)) {
    stop("iter must be a positive integer scalar.")
  }
  
  # Split one string into individual UTF-8 characters.
  split_chars <- function(x) {
    strsplit(enc2utf8(x), "", perl = TRUE)[[1]]
  }
  
  #check no empty strings
  if(any(nchar(S)==0)){
    stop("Sequences must not contain empty strings.")
  }
  
  # Convert each input string into a character vector.
  seq_chars <- lapply(S, split_chars)
  n <- length(seq_chars)
  lens <- vapply(seq_chars, length, integer(1))
  
  # Every sequence must be at least as long as the motif width.
  if (any(lens < w)) {
    stop("All sequences must have length >= w.")
  }
  
  # Build the alphabet from all characters observed across all sequences.
  alphabet <- sort(unique(unlist(seq_chars, use.names = FALSE)))
  a <- length(alphabet)
  char_index <- setNames(seq_len(a), alphabet)
  
  # Small pseudocount avoids zero probabilities in the profile.
  pseudocount <- 1
  
  # For each sequence, enumerate all valid start positions for a width-w window.
  valid_starts <- lapply(lens, function(L) seq_len(L - w + 1L))
  
  # Initialize each sequence with a random motif/window start.
  starts <- vapply(valid_starts, function(v) sample(v, 1L), integer(1))
  
  # Build a position-specific profile matrix from all windows except any excluded
  # sequence(s). Rows correspond to motif positions, columns to alphabet letters.
  build_profile <- function(starts, exclude = integer(0)) {
    counts <- matrix(
      pseudocount,
      nrow = w,
      ncol = a,
      dimnames = list(NULL, alphabet)
    )
    
    # Keep all sequences except those excluded for leave-one-out Gibbs updates.
    keep <- setdiff(seq_len(n), exclude)
    for (i in keep) {
      # Extract the current width-w window for sequence i.
      win <- seq_chars[[i]][starts[i]:(starts[i] + w - 1L)]
      cols <- unname(char_index[win])
      
      # Increment counts position by position.
      for (pos in seq_len(w)) {
        counts[pos, cols[pos]] <- counts[pos, cols[pos]] + 1
      }
    }
    
    # Convert counts to probabilities row-wise.
    counts / rowSums(counts)
  }
  
  # Compute the log-probability of choosing a given start in sequence i under
  # the provided profile.
  window_log_prob <- function(i, start, profile) {
    win <- seq_chars[[i]][start:(start + w - 1L)]
    cols <- unname(char_index[win])
    idx <- cbind(seq_len(w), cols)
    sum(log(profile[idx]))
  }
  
  # Gibbs update for one sequence:
  # 1) Build a profile from all other sequences.
  # 2) Score every candidate start.
  # 3) Sample a new start proportional to its probability.
  sample_start <- function(i, starts) {
    profile <- build_profile(starts, exclude = i)
    candidates <- valid_starts[[i]]
    logp <- vapply(candidates, function(s) window_log_prob(i, s, profile), numeric(1))
    
    # Stabilize exponentiation by subtracting the maximum log-probability.
    m <- max(logp)
    p <- exp(logp - m)
    p <- p / sum(p)
    sample(candidates, size = 1L, prob = p)
  }
  
  # Score a complete configuration of start positions by summing the log-
  # probabilities of all selected windows under the induced profile.
  alignment_score <- function(starts) {
    profile <- build_profile(starts)
    total <- 0
    for (i in seq_len(n)) {
      total <- total + window_log_prob(i, starts[i], profile)
    }
    total
  }
  
  # Track the best configuration observed during sampling.
  best_starts <- starts
  best_score <- alignment_score(starts)
  
  # Run Gibbs sweeps. In each sweep, update sequences in random order.
  for (it in seq_len(iter)) {
    for (i in sample.int(n)) {
      starts[i] <- sample_start(i, starts)
    }
    
    # Keep the highest-scoring alignment encountered.
    current_score <- alignment_score(starts)
    if (current_score > best_score) {
      best_score <- current_score
      best_starts <- starts
    }
  }
  
  # Reconstruct the full alignment while preserving unaligned ends.
  #
  # Idea:
  # - Each sequence keeps all of its original characters.
  # - We shift rows left/right so the sampled motif starts occur in the same
  #   alignment column.
  # - Any empty cells introduced by shifting are filled with '-'.
  anchor_col <- max(best_starts)
  left_pad <- anchor_col - best_starts
  row_widths <- left_pad + lens
  total_cols <- max(row_widths)
  
  # Initialize the rectangular output matrix with gap characters.
  out <- matrix("-", nrow = n, ncol = total_cols)
  rownames(out) <- names(S)
  
  # Write each original sequence into its shifted position.
  for (i in seq_len(n)) {
    cols <- (left_pad[i] + 1L):(left_pad[i] + lens[i])
    out[i, cols] <- seq_chars[[i]]
  }
  
  # Store useful metadata on the returned alignment.
  attr(out, "starts") <- best_starts
  attr(out, "motif_width") <- w
  attr(out, "score") <- best_score
  
  return(out)
}

#examples
S = c("XABBA", "XXABBA", "ABBAX")
gibbs_aligner_global(S, w = 4, iter = 10)

data = c("ABBACC", "BBABBA","CABBA")
output = gibbs_aligner_global(data, w = 4, iter = 10)

data = c("AAABBACC", "CCBBABBA","CABBA")
output = gibbs_aligner_global(data, w = 4, iter = 500)

data = c("ABBACCAA", "GGBBABBAG","CABBA", "GCAKCAABBTX")
output = gibbs_aligner_global(data, w = 4, iter = 10)
