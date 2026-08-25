library(parallel)
library(tibble)
library(readr)


## ============================================================
## PROGRESSIVE ALIGNMENT -- PARALLEL ANALYSIS
## ============================================================


## ------------------------------------------------------------
## 1. Tuning parameters
## ------------------------------------------------------------

match <- c(1, 2)
mismatch <- c(-1, -2)
gap <- c(-1, -2)

pa_tune_params <- expand.grid(
  match = match,
  mismatch = mismatch,
  gap = gap
)


## ------------------------------------------------------------
## 2. Function to process one lineage-pair
## ------------------------------------------------------------

run_pair <- function(i) {
  
  pair_dat <- paired_data[[i]]
  
  # Get Bird IDs for labelling alignment
  IDs <- pair_dat$Bird.ID
  bIDs <- paste0(IDs,"_",ave(IDs, IDs, FUN = seq_along))
  
  # Get note sequences
  bird_songsseqs <- pair_dat$note.seq
  
  # Store results for each parameter combination
  score_stats <- vector(
    "list",
    nrow(pa_tune_params)
  )
  
  for(k in seq_len(nrow(pa_tune_params))) {
    
    params <- pa_tune_params[k, ]
    
    # Fit progressive alignment
    dynam_model <- progressive_align(
      S = bird_songsseqs,
      match = params$match,
      mismatch = params$mismatch,
      gap = params$gap
    )
    
    # Reorder alignment matrix
    A <- og_order(
      align_mat = dynam_model,
      song_seqs = bird_songsseqs
    )
    
    # Add sequence names
    rownames(A) <- bIDs
    
    # Compute minimum entropy score
    line_scores <- min_entropy(A)
    align_len <- ncol(A)
    
    score_stats[[k]] <- tibble::tibble(
      score = line_scores,
      length = align_len,
      norm_score = line_scores / align_len
    )
  }
  
  # Gather score statistics
  scores <- do.call(
    rbind,
    score_stats
  )
  
  # Pair name
  pair_name <- names(paired_data)[i]
  
  # Parameter information
  res <- tibble::tibble(
    pa_tune_params,
    pair = pair_name
  )
  
  # Return results
  return(cbind(res, scores))
}


## ------------------------------------------------------------
## 3. Start cluster
## ------------------------------------------------------------

# 8 workers = 8 simultaneous lineage-pair jobs.
#
# If your machine has fewer than ~12 logical CPU cores,
# consider using 6 instead.

cl <- makeCluster(8)


## ------------------------------------------------------------
## 4. Load custom functions onto every worker
## ------------------------------------------------------------

clusterEvalQ(cl, {
  
  library(tibble)
  
  # progressive_align() sources:
  source("./functions/needleman.R")
  source("./functions/dynamic_program/compute_distance_matrix.R")
  source("./alignment_scripts/dynami_program/progressive_align.R")
  source("./functions/get_column_stats.R")
  source("./functions/min_entropy.R")
})

#checks

clusterEvalQ(cl, {
  c(
    tibble = requireNamespace("tibble", quietly = TRUE),
    progressive_align = exists("progressive_align"),
    needleman = exists("needleman"),
    min_entropy = exists("min_entropy")
  )
})

clusterEvalQ(cl, {
  exists("tibble")
})


## ------------------------------------------------------------
## 5. Export data, parameters and local functions
## ------------------------------------------------------------

clusterExport(
  cl,
  c(
    "paired_data",
    "pa_tune_params",
    "run_pair",
    "og_order"
  )
)


## ------------------------------------------------------------
## 6. TEST -- ONLY TWO PAIRS
## ------------------------------------------------------------

cat("Running test on pairs 1 and 2...\n")

test_results <- parLapply(
  cl,
  1:2,
  run_pair
)

cat("Test completed.\n")

#8 workers check
system.time({
  test_results <- parLapply(
    cl,
    1:8,
    run_pair
  )
})


## ------------------------------------------------------------
## 7. Inspect test output
## ------------------------------------------------------------

print(test_results[[1]])
print(test_results[[2]])


## ------------------------------------------------------------
## 8. FULL RUN -- ALL PAIRS
## ------------------------------------------------------------

cat(
  "Starting full progressive alignment run on",
  length(paired_data),
  "lineage pairs...\n"
)

start_time <- Sys.time()

dynam_scores <- parLapply(
  cl,
  seq_along(paired_data),
  run_pair
)

end_time <- Sys.time()

cat("Full run completed.\n")
cat("Started:", start_time, "\n")
cat("Finished:", end_time, "\n")
cat("Elapsed:", end_time - start_time, "\n")


## ------------------------------------------------------------
## 9. Stop cluster
## ------------------------------------------------------------

stopCluster(cl)


## ------------------------------------------------------------
## 10. Combine results
## ------------------------------------------------------------

final_dp_res <- do.call(
  rbind,
  dynam_scores
)


## ------------------------------------------------------------
## 11. Save results
## ------------------------------------------------------------

readr::write_csv(
  final_dp_res,
  file = "./results/msa_scores/inter_dp_scores.csv"
)

cat("Results written to:\n")
cat("./results/msa_scores/inter_dp_scores.csv\n")