#helper functions for gibbs_aligner_global

#do not export

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

#tests

# Compute the log-probability of choosing a given start in sequence i under
# the provided profile.
window_log_prob <- function(i, start, profile) {
  win <- seq_chars[[i]][start:(start + w - 1L)]
  cols <- unname(char_index[win])
  idx <- cbind(seq_len(w), cols)
  sum(log(profile[idx]))
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

#tests----

library(testthat)

library(testthat)

test_that("build_profile creates valid profile matrix", {
  
  # recreate variables from gibbs_aligner_global()
  S <- c("ABCX", "ABCY", "ABCZ")
  
  seq_chars <- lapply(S, function(x) strsplit(x, "")[[1]])
  
  n <- length(seq_chars)
  w <- 3
  
  alphabet <- sort(unique(unlist(seq_chars)))
  a <- length(alphabet)
  
  char_index <- setNames(seq_along(alphabet), alphabet)
  
  pseudocount <- 1
  
  
  # Create environment containing required variables
  env <- new.env()
  
  list2env(
    list(
      seq_chars = seq_chars,
      n = n,
      w = w,
      a = a,
      alphabet = alphabet,
      char_index = char_index,
      pseudocount = pseudocount
    ),
    envir = env
  )
  
  # Put function into that environment
  environment(build_profile) <- env
  
  
  profile <- build_profile(
    starts = c(1,1,1)
  )
  
  
  expect_equal(
    dim(profile),
    c(3,length(alphabet))
  )
  
  expect_equal(
    rowSums(profile),
    rep(1,3)
  )
  
  consensus <- alphabet[max.col(profile)]
  
  #check final motif alignment
  expect_equal(
    consensus,
    c("A","B","C")
  )
})

test_that("window_log_prob prefers conserved windows", {
  
  S <- c(
    "XXABC",
    "YYABC",
    "ZZABC"
  )
  
  seq_chars <- lapply(S, function(x) strsplit(x,"")[[1]])
  
  n <- length(seq_chars)
  w <- 3
  
  alphabet <- sort(unique(unlist(seq_chars)))
  a <- length(alphabet)
  
  char_index <- setNames(seq_along(alphabet), alphabet)
  
  pseudocount <- 1
  
  
  env <- new.env()
  
  list2env(
    list(
      seq_chars = seq_chars,
      n = n,
      w = w,
      a = a,
      alphabet = alphabet,
      char_index = char_index,
      pseudocount = pseudocount
    ),
    envir = env
  )
  
  
  environment(build_profile) <- env
  environment(window_log_prob) <- env
  
  
  profile <- build_profile(
    starts=c(3,3,3),
    exclude=3
  )
  
  
  good <- window_log_prob(
    i=3,
    start=3,
    profile=profile
  )
  
  bad <- window_log_prob(
    i=3,
    start=1,
    profile=profile
  )
  
  
  expect_gt(good,bad)
  
})

test_that("alignment_score favours conserved alignment", {
  
  S <- c(
    "XXABC",
    "YYABC",
    "ZZABC"
  )
  
  seq_chars <- lapply(S, function(x) strsplit(x,"")[[1]])
  
  n <- length(seq_chars)
  w <- 3
  
  alphabet <- sort(unique(unlist(seq_chars)))
  a <- length(alphabet)
  
  char_index <- setNames(seq_along(alphabet), alphabet)
  
  pseudocount <- 1
  
  
  env <- new.env()
  
  list2env(
    list(
      seq_chars = seq_chars,
      n = n,
      w = w,
      a = a,
      alphabet = alphabet,
      char_index = char_index,
      pseudocount = pseudocount
    ),
    envir = env
  )
  
  
  environment(build_profile) <- env
  environment(window_log_prob) <- env
  environment(alignment_score) <- env
  
  
  good <- alignment_score(
    c(3,3,3)
  )
  
  bad <- alignment_score(
    c(1,1,1)
  )
  
  
  expect_gt(
    good,
    bad
  )
  
})
