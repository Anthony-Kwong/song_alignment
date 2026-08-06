#' Create a balanced inter-lineage sample.
#'
#' Randomly samples an equal number of observations from two lineages. The
#' larger lineage is downsampled to match the size of the smaller lineage.
#'
#' @param df Data frame containing lineage information.
#' @param line1 First lineage.
#' @param line2 Second lineage.
#' @param line_col Name of the lineage column.
#'
#' @return A data frame containing a balanced sample from the two lineages.
sample_lineage_pair <- function(df,
                                line1,
                                line2,
                                line_col = "Line") {
  
  df1 <- df[df[[line_col]] == line1, ]
  df2 <- df[df[[line_col]] == line2, ]
  
  n_sample <- min(nrow(df1), nrow(df2))
  
  df1 <- df1[sample(seq_len(nrow(df1)), n_sample), , drop = FALSE]
  df2 <- df2[sample(seq_len(nrow(df2)), n_sample), , drop = FALSE]
  
  rbind(df1, df2)
}


#' Generate balanced samples for every lineage pair.
#'
#' Creates balanced inter-lineage datasets for all pairwise lineage
#' combinations.
#'
#' @param df Data frame containing lineage information.
#' @param line_col Name of the lineage column.
#'
#' @return A named list of balanced data frames.
all_lineage_pairs <- function(df,
                              line_col = "Line") {
  
  lines <- unique(df[[line_col]])
  pairs <- combn(lines, 2, simplify = FALSE)
  out <- lapply(
    pairs,
    function(x) {
      sample_lineage_pair(
        df,
        line1 = x[1],
        line2 = x[2],
        line_col = line_col
      )
    }
  )
  names(out) <- vapply(
    pairs,
    function(x) paste(x, collapse = "_vs_"),
    character(1)
  )
  return(out)
}

#tests 

library(testthat)
test_that("sample_lineage_pair returns equal numbers from each lineage", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 5),
             rep("B", 3),
             rep("C", 4)),
    Song = seq_len(12)
  )
  
  x <- sample_lineage_pair(df, "A", "B")
  
  expect_equal(nrow(x), 6)
 
  expect_equal(
    as.numeric(table(x$Line)),
    c(3, 3)
  )
})


test_that("sample_lineage_pair never samples more than available", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 2),
             rep("B", 7)),
    Song = seq_len(9)
  )
  
  x <- sample_lineage_pair(df, "A", "B")
  
  expect_equal(nrow(x), 4)
  
  expect_equal(
    as.numeric(table(x$Line)),
    c(2, 2)
  )
  
})


test_that("sample_lineage_pair only returns requested lineages", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 5),
             rep("B", 4),
             rep("C", 6)),
    Song = seq_len(15)
  )
  
  x <- sample_lineage_pair(df, "A", "C")
  
  expect_true(
    all(x$Line %in% c("A", "C"))
  )
  
})


test_that("all_lineage_pairs returns every pair", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 5),
             rep("B", 3),
             rep("C", 4)),
    Song = seq_len(12)
  )
  
  x <- all_lineage_pairs(df)
  
  expect_equal(length(x), 3)
  
  expect_equal(
    sort(names(x)),
    sort(c(
      "A_vs_B",
      "A_vs_C",
      "B_vs_C"
    ))
  )
  
})


test_that("all returned lineage pairs are balanced", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 6),
             rep("B", 3),
             rep("C", 5)),
    Song = seq_len(14)
  )
  
  x <- all_lineage_pairs(df)
  
  for(i in seq_along(x)){
    
    counts <- table(x[[i]]$Line)
    
    expect_equal(length(counts), 2)
    
    expect_equal(
      as.numeric(counts[1]),
      as.numeric(counts[2])
    )
    
  }
  
})


test_that("sampled rows always come from original dataframe", {
  
  set.seed(1)
  
  df <- data.frame(
    Line = c(rep("A", 5),
             rep("B", 4)),
    Song = seq_len(9)
  )
  
  x <- sample_lineage_pair(df, "A", "B")
  
  expect_true(
    all(x$Song %in% df$Song)
  )
  
})

