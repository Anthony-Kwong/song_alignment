#produce results tables and figures for alignment comparison

#library
library(dplyr)
library(ggplot2)


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


# Existing files have one contiguous block per paired_data[[i]]:
# PHMM: 9 rows/run; Gibbs: 6 rows/run; progressive: 8 rows/run.
# The same paired_data order is assumed across all three methods.

add_run_id <- function(x, n_params, method) {
  n_pairs <- n_distinct(x$pair)
  rows_per_run <- n_pairs * n_params
  
  x %>%
    mutate(
      method = method,
      run_id = (row_number() - 1L) %/% rows_per_run + 1L
    )
}

#adding run ids ----
inter_phmm <- add_run_id(inter_phmm, 9, "PHMM")
inter_gibbs <- add_run_id(inter_gibbs, 6, "Gibbs")
inter_progressive <- add_run_id(inter_pa, 8, "Progressive")

#checks ----


# Check that every method recovered the same number and order of runs.
run_check <- bind_rows(
  inter_phmm %>% distinct(method, run_id, pair),
  inter_gibbs %>% distinct(method, run_id, pair),
  inter_progressive %>% distinct(method, run_id, pair)
) %>%
  arrange(run_id, method)

print(run_check %>% count(method)) #check agreement

run_check %>%
  count(method, run_id, name = "n_pairs")

run_check %>%
  count(method, run_id, name = "n_pairs") %>%
  count(n_pairs, name = "n_method_runs")

# Check how many sampling replicates exist per lineage pair.
# If the result is 5 rather than 10, the saved data contain 5 runs/pair.
replicate_check <- run_check %>%
  group_by(method, pair) %>%
  summarise(
    n_sampling_replicates = n_distinct(run_id),
    .groups = "drop"
  )

print(replicate_check %>% count(method, n_sampling_replicates))

#results pooled across input params, per method, per pair

#ensures cohens D is from 10 intra scores, vs 10 choose 2 inter scores

phmm_run_scores <- inter_phmm %>%
  group_by(pair, run_id) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

#get mean scores per pair
phmm_pair_scores <- phmm_run_scores %>%
  group_by(pair) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

gibbs_run_scores <-inter_gibbs %>%
  group_by(pair, run_id) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

gibbs_pair_scores <- gibbs_run_scores %>%
  group_by(pair) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

progressive_run_scores <- inter_progressive %>%
  group_by(pair, run_id) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

progressive_pair_scores <- progressive_run_scores %>%
  group_by(pair) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

#average results across params, per lineage for intra-lineage results

# Average parameterisations within each lineage
intra_gibbs_final <- intra_gibbs %>%
  group_by(line) %>%                 # replace pair with your lineage ID column if needed
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

intra_pa_final <- intra_pa %>%
  group_by(line) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )

intra_phmm_final <- intra_phmm %>%
  group_by(comp) %>%
  summarise(
    norm_score = mean(norm_score, na.rm = TRUE),
    .groups = "drop"
  )


#main result table (global, pooled params) ----

library(dplyr)
library(tibble)

make_method_row <- function(method, intra, inter) {

  tibble(
    method = method,
    n_intra = length(intra),
    mean_intra = mean(intra, na.rm = TRUE),
    sd_intra = sd(intra, na.rm = TRUE),
    n_inter = length(inter),
    mean_inter = mean(inter, na.rm = TRUE),
    sd_inter = sd(inter, na.rm = TRUE),
    mean_difference = mean(inter, na.rm = TRUE) - mean(intra, na.rm = TRUE),
    cohen_d = cohens_d(intra = intra,inter = inter)
  )
}

main_results <- bind_rows(
  make_method_row(
    "Gibbs",
    intra = intra_gibbs_final$norm_score,
    inter = gibbs_pair_scores$norm_score
  ),
  make_method_row(
    "Progressive",
    intra = intra_pa_final$norm_score,
    inter = progressive_pair_scores$norm_score
  ),
  make_method_row(
    "PHMM",
    intra = intra_phmm_final$norm_score,
    inter = phmm_pair_scores$norm_score
  )
) %>%
  #round to 3dp
  mutate(
    across(
      c(mean_intra, sd_intra, mean_inter, sd_inter,
        mean_difference, cohen_d),
      ~ round(.x, 3)
    )
  )

print(main_results)

#select columns for final table
main_tab = main_results %>%
  dplyr::select(method, mean_intra, sd_intra, mean_inter, sd_inter, cohen_d)

main_tab <- main_tab |>
  rename(
    Method = method,
    `Intra-lineage mean` = mean_intra,
    `Intra-lineage SD` = sd_intra,
    `Inter-lineage mean` = mean_inter,
    `Inter-lineage SD` = sd_inter,
    `Cohen's $d$` = cohen_d
  )

print(xtable::xtable(main_tab), file = "./results/msa_scores/global_comp.tex")

#boxplot of global results ----

plot_dat <- bind_rows(
  tibble(
    method = "Gibbs",
    type = "Intra-lineage",
    score = intra_gibbs_final$norm_score
  ),
  tibble(
    method = "Gibbs",
    type = "Inter-lineage",
    score = gibbs_pair_scores$norm_score
  ),
  tibble(
    method = "Progressive",
    type = "Intra-lineage",
    score = intra_pa_final$norm_score
  ),
  tibble(
    method = "Progressive",
    type = "Inter-lineage",
    score = progressive_pair_scores$norm_score
  ),
  tibble(
    method = "PHMM",
    type = "Intra-lineage",
    score = intra_phmm_final$norm_score
  ),
  tibble(
    method = "PHMM",
    type = "Inter-lineage",
    score = phmm_pair_scores$norm_score
  )
)

plot_dat$method <- factor(
  plot_dat$method,
  levels = c("Gibbs", "Progressive", "PHMM")
)

plot_dat$type <- factor(
  plot_dat$type,
  levels = c("Intra-lineage", "Inter-lineage")
)

p <- ggplot(plot_dat, aes(x = method, y = score, fill = type)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(group = type),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    alpha = 0.35,
    size = 1.5
  ) +
  labs(
    x = "Alignment method",
    y = "Normalised entropy",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(p, file ="./results/msa_scores/global_boxplot.png")

#sensitivity analysis plots ----

# Average the 5 repeated runs for each lineage pair
inter_phmm_summary <- inter_phmm %>%
  group_by(lambda, max_scale, pair) %>%
  summarise(
    norm_score = mean(norm_score),
    .groups = "drop"
  )

# Calculate Cohen's d for each PHMM parameter combination
phmm_sensitivity <- intra_phmm %>%
  
  # Group intra-lineage scores by the two PHMM input parameters
  group_by(lambda, max_scale) %>%
  
  # Calculate the effect size for each parameter combination
  group_modify(~ {
    
    # Extract the 10 intra-lineage scores for this parameter combination
    intra <- .x$norm_score
    
    # Extract the 45 inter-lineage scores for the same parameter combination
    inter <- inter_phmm_summary %>%
      filter(
        lambda == .y$lambda,
        max_scale == .y$max_scale
      ) %>%
      pull(norm_score)
    
    # Calculate Cohen's d, with positive values indicating
    # higher entropy between lineages than within lineages
    tibble(
      cohen_d = cohens_d(
        intra = intra,
        inter = inter
      )
    )
    
  }) %>%
  
  # Remove grouping from the final sensitivity results
  ungroup()

#check

intra_test <- intra_phmm %>%
  filter(lambda == 0, max_scale == 1) %>%
  pull(norm_score)

inter_test <- inter_phmm_summary %>%
  filter(lambda == 0, max_scale == 1) %>%
  pull(norm_score)

length(intra_test) # 10
length(inter_test) # 45, 10 C 2
#check agreement
cohens_d(intra = intra_test, inter = inter_test)

# Average repeated runs for each lineage pair and parameterisation
inter_progressive_summary <- inter_progressive %>%
  group_by(match, mismatch, gap, pair) %>%
  summarise(
    norm_score = mean(norm_score),
    .groups = "drop"
  )

# Calculate Cohen's d for each parameter combination
progressive_sensitivity <- intra_pa %>%
  group_by(match, mismatch, gap) %>%
  group_modify(~ {
    
    # Extract intra-lineage scores for this parameter combination
    intra <- .x$norm_score
    
    # Extract inter-lineage scores for the same parameter combination
    inter <- inter_progressive_summary %>%
      filter(
        match == .y$match,
        mismatch == .y$mismatch,
        gap == .y$gap
      ) %>%
      pull(norm_score)
    
    # Calculate Cohen's d
    tibble(
      cohen_d = cohens_d(
        intra = intra,
        inter = inter
      )
    )
    
  }) %>%
  ungroup()

progressive_sensitivity

#check

intra_test <- intra_pa %>%
  filter(match == 1, mismatch == -2, gap == -1) %>%
  pull(norm_score)

inter_test <- inter_progressive_summary %>%
  filter(match == 1, mismatch == -2, gap == -1) %>%
  pull(norm_score)

length(intra_test) # 10
length(inter_test) # 45, 10 C 2
#check agreement
cohens_d(intra = intra_test, inter = inter_test)

# Average the repeated runs for each lineage pair and wminus value
inter_gibbs_summary <- inter_gibbs %>%
  group_by(wminus, pair) %>%
  summarise(
    norm_score = mean(norm_score),
    .groups = "drop"
  )

# Calculate Cohen's d for each wminus value
gibbs_sensitivity <- intra_gibbs %>%
  group_by(wminus) %>%
  group_modify(~ {
    
    # Extract intra-lineage scores for this parameter value
    intra <- .x$norm_score
    
    # Extract inter-lineage scores for the same parameter value
    inter <- inter_gibbs_summary %>%
      filter(wminus == .y$wminus) %>%
      pull(norm_score)
    
    # Calculate Cohen's d
    tibble(
      cohen_d = cohens_d(
        intra = intra,
        inter = inter
      )
    )
    
  }) %>%
  ungroup()

#check

intra_test <- intra_gibbs %>%
  filter(wminus == 0) %>%
  pull(norm_score)

inter_test <- inter_gibbs_summary %>%
  filter(wminus == 0) %>%
  pull(norm_score)

length(intra_test) # 10
length(inter_test) # 45, 10 C 2
#check agreement
cohens_d(intra = intra_test, inter = inter_test)

gibbs_sensitivity

#gibbs sensitivity plot ----
gibbs_plot <- ggplot(
  gibbs_sensitivity,
  aes(x = wminus, y = cohen_d)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(limits = c(0.55, 0.80)) +
  labs(
    x = expression(w[minus]),
    y = "Cohen's d"
  ) +
  theme_classic(base_size = 12)

#phmm sensitivity plot ----
phmm_plot <- ggplot(
  phmm_sensitivity,
  aes(
    x = factor(max_scale),
    y = factor(lambda),
    fill = cohen_d
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", cohen_d)),
    size = 4
  ) +
  scale_fill_viridis_c(
    limits = c(0.55, 0.80)
  ) +
  labs(
    x = expression(max~scale),
    y = expression(lambda),
    fill = "Cohen's d"
  ) +
  theme_classic(base_size = 12)

#progressive sensitivity plot
progressive_plot <- ggplot(
  progressive_sensitivity,
  aes(
    x = factor(mismatch),
    y = factor(match),
    fill = cohen_d
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", cohen_d)),
    size = 4
  ) +
  scale_fill_viridis_c(
    limits = c(0.55, 0.80)
  ) +
  facet_wrap(
    ~ gap,
    labeller = label_bquote(gap == .(gap))
  ) +
  labs(
    x = "Mismatch score",
    y = "Match score",
    fill = "Cohen's d"
  ) +
  theme_classic(base_size = 12)

library(ggpubr)

sensitivity_figure <- ggarrange(
  gibbs_plot,
  progressive_plot,
  phmm_plot,
  ncol = 1,
  labels = c("A", "B", "C"),
  align = "v"
)

ggsave(
  filename = "./results/msa_scores/sensitivity_figure.pdf",
  plot = sensitivity_figure,
  width = 6.5,
  height = 8.5,
  units = "in",
  device = cairo_pdf
)

#sensitivity table

sensitivity_summary <- tibble(
  method = c("Gibbs", "Progressive", "PHMM"),
  
  overall_d = c(
    main_results$cohen_d[main_results$method == "Gibbs"],
    main_results$cohen_d[main_results$method == "Progressive"],
    main_results$cohen_d[main_results$method == "PHMM"]
  ),
  
  min_d = c(
    min(gibbs_sensitivity$cohen_d),
    min(progressive_sensitivity$cohen_d),
    min(phmm_sensitivity$cohen_d)
  ),
  
  max_d = c(
    max(gibbs_sensitivity$cohen_d),
    max(progressive_sensitivity$cohen_d),
    max(phmm_sensitivity$cohen_d)
  )
) %>%
  mutate(
    overall_d = sprintf("%.2f", overall_d),
    sensitivity_range = sprintf("%.2f--%.2f", min_d, max_d)
  ) %>%
  select(
    Method = method,
    `Overall Cohen's d` = overall_d,
    `Range across parameterisations` = sensitivity_range
  )

sensitivity_summary

print(
  xtable::xtable(
    sensitivity_summary,
    caption = paste(
      "Overall Cohen's d and range of Cohen's d across the",
      "parameterisations tested for each alignment method."
    ),
    label = "tab:sensitivity_summary"
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)
