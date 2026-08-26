#produce results tables and figures for alignment comparison

#library
library(dplyr)

#load inter-lineage results
inter_phmm = readr::read_csv("./results/msa_scores/inter_phmm_scores.csv")
inter_gibbs = readr::read_csv("./results/msa_scores/inter_gibbs_scores.csv")
inter_pa = readr::read_csv("./results/msa_scores/inter_dp_scores.csv")

#load intra-lineage results
intra_phmm = readr::read_csv("./results/msa_scores/phmm_scores.csv")
intra_gibbs = readr::read_csv("./results/msa_scores/gibbs_scores.csv")
intra_pa = readr::read_csv("./results/msa_scores/dp_scores.csv")

# Cohen's d: positive = inter-lineage scores are higher

#compute cohen'D using intra- and inter-lineage scores
#intra: numeric vector of intra scores
#inter: numeric vector of inter scores
#output cohen's D score

cohens_d <- function(intra, inter) {
  n1 <- length(intra)
  n2 <- length(inter)
  
  pooled_sd <- sqrt(
    ((n1 - 1) * var(intra) +
       (n2 - 1) * var(inter)) /
      (n1 + n2 - 2)
  )
  
  res = (mean(inter) - mean(intra)) / pooled_sd
  return(res)
}

#computing scores per input params ----

gibbs_results <- bind_rows(
  intra_gibbs %>% mutate(type = "Intra"),
  inter_gibbs %>% mutate(type = "Inter")
) %>%
  group_by(wminus) %>%
  summarise(
    n_intra = sum(type == "Intra"),
    n_inter = sum(type == "Inter"),
    
    mean_intra = mean(norm_score[type == "Intra"]),
    sd_intra = sd(norm_score[type == "Intra"]),
    
    mean_inter = mean(norm_score[type == "Inter"]),
    sd_inter = sd(norm_score[type == "Inter"]),
    
    cohen_d = cohens_d(
      norm_score[type == "Intra"],
      norm_score[type == "Inter"]
    ),
    
    .groups = "drop"
  )

phmm_results <- bind_rows(
  intra_phmm %>% mutate(type = "Intra"),
  inter_phmm %>% mutate(type = "Inter")
) %>%
  group_by(lambda, max_scale) %>%
  summarise(
    n_intra = sum(type == "Intra"),
    n_inter = sum(type == "Inter"),
    
    mean_intra = mean(norm_score[type == "Intra"]),
    sd_intra = sd(norm_score[type == "Intra"]),
    
    mean_inter = mean(norm_score[type == "Inter"]),
    sd_inter = sd(norm_score[type == "Inter"]),
    
    cohen_d = cohens_d(
      norm_score[type == "Intra"],
      norm_score[type == "Inter"]
    ),
    
    .groups = "drop"
  )

progressive_results <- bind_rows(
  intra_pa %>% mutate(type = "Intra"),
  inter_pa %>% mutate(type = "Inter")
) %>%
  group_by(match, mismatch, gap) %>%
  summarise(
    n_intra = sum(type == "Intra"),
    n_inter = sum(type == "Inter"),
    
    mean_intra = mean(norm_score[type == "Intra"]),
    sd_intra = sd(norm_score[type == "Intra"]),
    
    mean_inter = mean(norm_score[type == "Inter"]),
    sd_inter = sd(norm_score[type == "Inter"]),
    
    cohen_d = cohens_d(
      norm_score[type == "Intra"],
      norm_score[type == "Inter"]
    ),
    
    .groups = "drop"
  )

