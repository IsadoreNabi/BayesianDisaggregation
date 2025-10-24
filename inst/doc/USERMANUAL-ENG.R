## ----setup, include=TRUE------------------------------------------------------
# Global chunk defaults
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 9, fig.height = 6
)

# Libraries
suppressPackageStartupMessages({
  library(BayesianDisaggregation)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(openxlsx)
})

# Logging verbosity from the package
log_enable("INFO")
set.seed(2024)

## ----L_from_P, include=TRUE---------------------------------------------------
# Example call (internals are in the package):
# L <- compute_L_from_P(P) 

## ----spread-L, include=TRUE---------------------------------------------------
# Example call:
# LT <- spread_likelihood(L, T_periods = nrow(P), pattern = "recent")

## ----posteriors, include=TRUE-------------------------------------------------
# posterior_weighted(P, LT, lambda = 0.7)
# posterior_multiplicative(P, LT)
# posterior_dirichlet(P, LT, gamma = 0.1)
# posterior_adaptive(P, LT)

## ----metrics-fns, include=TRUE------------------------------------------------
# coherence_score(P, W, L, mult = 3.0, const = 0.5)
# numerical_stability_exp(W, a = 1000, b = 10)
# temporal_stability(W, kappa = 50)
# stability_composite(W, a = 1000, b = 10, kappa = 50)

## ----interp-fn, include=TRUE--------------------------------------------------
# interpretability_score(P, W, use_q90 = TRUE)

## ----api, include=TRUE--------------------------------------------------------
# Example signature (see Section 8 for real data):
# bayesian_disaggregate(path_cpi, path_weights,
#   method = c("weighted","multiplicative","dirichlet","adaptive"),
#   lambda = 0.7, gamma = 0.1,
#   coh_mult = 3.0, coh_const = 0.5,
#   stab_a = 1000, stab_b = 10, stab_kappa = 50,
#   likelihood_pattern = "recent")

## ----demo, include=TRUE-------------------------------------------------------
# Synthetic prior matrix (rows on simplex)
T <- 10; K <- 6
set.seed(123)
P <- matrix(rexp(T*K), nrow = T)
P <- P / rowSums(P)

# Likelihood vector from P (PCA/SVD; robust with fallback)
L  <- compute_L_from_P(P)

# Spread over time with "recent" pattern
LT <- spread_likelihood(L, T_periods = T, pattern = "recent")

# Try a couple of posteriors
W_weighted <- posterior_weighted(P, LT, lambda = 0.7)
W_adaptive <- posterior_adaptive(P, LT)

# Metrics for adaptive
coh  <- coherence_score(P, W_adaptive, L)
stab <- stability_composite(W_adaptive, a = 1000, b = 10, kappa = 50)
intr <- interpretability_score(P, W_adaptive)
eff  <- 0.65
comp <- 0.30*coh + 0.25*stab + 0.25*intr + 0.20*eff

data.frame(coherence = coh, stability = stab, interpretability = intr,
           efficiency = eff, composite = comp) %>% round(4)

## ----real-pipeline, eval=FALSE------------------------------------------------
# # === Create synthetic data for CRAN-compliant demo ===
# demo_dir <- tempdir()
# 
# # Create synthetic CPI data (mimicking your structure)
# set.seed(123)
# cpi_demo <- data.frame(
#   Year = 2000:2010,
#   CPI = cumsum(c(100, rnorm(10, 0.5, 2)))  # Random walk starting at 100
# )
# cpi_file <- file.path(demo_dir, "synthetic_cpi.xlsx")
# openxlsx::write.xlsx(cpi_demo, cpi_file)
# 
# # Create synthetic weights matrix (mimicking VAB weights structure)
# set.seed(456)
# years <- 2000:2010
# sectors <- c("Agriculture", "Manufacturing", "Services", "Construction", "Mining")
# 
# weights_demo <- data.frame(Year = years)
# for(sector in sectors) {
#   weights_demo[[sector]] <- runif(length(years), 0.05, 0.35)
# }
# # Normalize rows to sum to 1 (simplex constraint)
# weights_demo[, -1] <- weights_demo[, -1] / rowSums(weights_demo[, -1])
# weights_file <- file.path(demo_dir, "synthetic_weights.xlsx")
# openxlsx::write.xlsx(weights_demo, weights_file)
# 
# # Use synthetic data paths
# path_cpi <- cpi_file
# path_w <- weights_file
# out_dir <- demo_dir
# 
# cat("Using synthetic data for CRAN demo:\n")
# cat("CPI file:", path_cpi, "\n")
# cat("Weights file:", path_w, "\n")
# cat("Output directory:", out_dir, "\n")
# 
# # --- Base run (robust defaults) ---
# base_res <- bayesian_disaggregate(
#   path_cpi           = path_cpi,
#   path_weights       = path_w,
#   method             = "adaptive",
#   lambda             = 0.7,   # recorded in metrics; not used by "adaptive"
#   gamma              = 0.1,
#   coh_mult           = 3.0,
#   coh_const          = 0.5,
#   stab_a             = 1000,
#   stab_b             = 10,
#   stab_kappa         = 60,
#   likelihood_pattern = "recent"
# )
# xlsx_base <- save_results(base_res, out_dir = file.path(out_dir, "base"))
# print(base_res$metrics)
# 
# # --- Minimal grid search for demo (reduced size) ---
# n_cores <- 1  # Use single core for CRAN compliance
# grid_df <- expand.grid(
#   method             = c("weighted", "adaptive"),  # Reduced methods
#   lambda             = c(0.5, 0.7),               # Reduced options
#   gamma              = 0.1,                       # Single option
#   coh_mult           = 3.0,                       # Single option
#   coh_const          = 0.5,                       # Single option
#   stab_a             = 1000,
#   stab_b             = 10,
#   stab_kappa         = 60,                        # Single option
#   likelihood_pattern = "recent",                  # Single option
#   KEEP.OUT.ATTRS     = FALSE,
#   stringsAsFactors   = FALSE
# )
# 
# grid_res <- run_grid_search(
#   path_cpi     = path_cpi,
#   path_weights = path_w,
#   grid_df      = grid_df,
#   n_cores      = n_cores
# )
# write.csv(grid_res, file.path(out_dir, "grid_results.csv"), row.names = FALSE)
# 
# best_row <- grid_res %>% arrange(desc(composite)) %>% slice(1)
# print("Best configuration from grid search:")
# print(best_row)
# 
# # --- Re-run the best configuration for clean export ---
# best_res <- bayesian_disaggregate(
#   path_cpi           = path_cpi,
#   path_weights       = path_w,
#   method             = best_row$method,
#   lambda             = if (!is.na(best_row$lambda)) best_row$lambda else 0.7,
#   gamma              = if (!is.na(best_row$gamma))  best_row$gamma  else 0.1,
#   coh_mult           = best_row$coh_mult,
#   coh_const          = best_row$coh_const,
#   stab_a             = best_row$stab_a,
#   stab_b             = best_row$stab_b,
#   stab_kappa         = best_row$stab_kappa,
#   likelihood_pattern = best_row$likelihood_pattern
# )
# xlsx_best <- save_results(best_res, out_dir = file.path(out_dir, "best"))
# 
# # --- One Excel with everything (including hyperparameters) ---
# sector_summary <- tibble(
#   Sector          = colnames(best_res$posterior)[-1],
#   prior_mean      = colMeans(as.matrix(best_res$prior[, -1])),
#   posterior_mean  = colMeans(as.matrix(best_res$posterior[, -1]))
# )
# 
# wb <- createWorkbook()
# addWorksheet(wb, "Hyperparameters"); writeData(wb, "Hyperparameters", best_row)
# addWorksheet(wb, "Metrics");         writeData(wb, "Metrics", best_res$metrics)
# addWorksheet(wb, "Prior_P");         writeData(wb, "Prior_P", best_res$prior)
# addWorksheet(wb, "Posterior_W");     writeData(wb, "Posterior_W", best_res$posterior)
# addWorksheet(wb, "Likelihood_t");    writeData(wb, "Likelihood_t", best_res$likelihood_t)
# addWorksheet(wb, "Likelihood_L");    writeData(wb, "Likelihood_L", best_res$likelihood)
# addWorksheet(wb, "Sector_Summary");  writeData(wb, "Sector_Summary", sector_summary)
# 
# for (sh in c("Hyperparameters","Metrics","Prior_P","Posterior_W",
#              "Likelihood_t","Likelihood_L","Sector_Summary")) {
#   freezePane(wb, sh, firstRow = TRUE)
#   addFilter(wb, sh, rows = 1, cols = 1:ncol(readWorkbook(wb, sh)))
#   setColWidths(wb, sh, cols = 1:200, widths = "auto")
# }
# 
# # --- Add sectoral CPI (aggregate times posterior weights) ---
# W_post <- best_res$posterior           # Year + sectors
# cpi_df <- read_cpi(path_cpi)           # Year, CPI
# sector_cpi <- dplyr::left_join(W_post, cpi_df, by = "Year") %>%
#   dplyr::mutate(dplyr::across(-c(Year, CPI), ~ .x * CPI))
# 
# # Quality check: sector sums vs CPI
# check_sum <- sector_cpi %>%
#   dplyr::mutate(row_sum = rowSums(dplyr::across(-c(Year, CPI))),
#                 diff    = CPI - row_sum)
# cat("Quality check (first 5 rows):\n")
# print(head(check_sum, 5))
# 
# addWorksheet(wb, "Sector_CPI")
# writeData(wb, "Sector_CPI", sector_cpi)
# freezePane(wb, "Sector_CPI", firstRow = TRUE)
# addFilter(wb, "Sector_CPI", rows = 1, cols = 1:ncol(sector_cpi))
# setColWidths(wb, "Sector_CPI", cols = 1:200, widths = "auto")
# 
# excel_onefile <- file.path(out_dir, "best", "Best_Full_Output_withSectorCPI.xlsx")
# saveWorkbook(wb, excel_onefile, overwrite = TRUE)
# cat("Full results saved to:", excel_onefile, "\n")
# 
# # --- Quick plots (saved as PNGs) ---
# dir_plots <- file.path(out_dir, "best", "plots")
# if (!dir.exists(dir_plots)) dir.create(dir_plots, recursive = TRUE)
# 
# W_long <- best_res$posterior %>%
#   pivot_longer(-Year, names_to = "Sector", values_to = "Weight")
# p_heat <- ggplot(W_long, aes(Year, Sector, fill = Weight)) +
#   geom_tile() + scale_fill_viridis_c() +
#   labs(title = "Posterior weights (W): heatmap", x = "Year", y = "Sector", fill = "Share") +
#   theme_minimal(base_size = 11) + theme(axis.text.y = element_text(size = 6))
# ggsave(file.path(dir_plots, "posterior_heatmap.png"), p_heat, width = 12, height = 9, dpi = 220)
# 
# top_sectors <- best_res$posterior %>%
#   summarise(across(-Year, mean)) %>%
#   pivot_longer(everything(), names_to = "Sector", values_to = "MeanShare") %>%
#   arrange(desc(MeanShare)) %>% slice(1:3) %>% pull(Sector)  # Reduced to top 3 for demo
# 
# p_lines <- best_res$posterior %>%
#   select(Year, all_of(top_sectors)) %>%
#   pivot_longer(-Year, names_to = "Sector", values_to = "Weight") %>%
#   ggplot(aes(Year, Weight, color = Sector)) +
#   geom_line(linewidth = 0.9) +
#   labs(title = "Top 3 sectors by average share (posterior W)", y = "Share", x = "Year") +
#   theme_minimal(base_size = 11)
# ggsave(file.path(dir_plots, "posterior_topSectors.png"), p_lines, width = 11, height = 6, dpi = 220)
# 
# cat("Demo completed successfully. All files written to temporary directory.\n")

## ----invariants, include=TRUE-------------------------------------------------
# Example: invariants on a fresh synthetic run
T <- 6; K <- 5
set.seed(7)
P <- matrix(rexp(T*K), nrow = T); P <- P / rowSums(P)
L <- compute_L_from_P(P)
LT <- spread_likelihood(L, T, "recent")
W  <- posterior_multiplicative(P, LT)

# Invariants
stopifnot(all(abs(rowSums(P)  - 1) < 1e-12))
stopifnot(all(abs(rowSums(LT) - 1) < 1e-12))
stopifnot(all(abs(rowSums(W)  - 1) < 1e-12))
c(
  coherence = coherence_score(P, W, L),
  stability = stability_composite(W),
  interpret = interpretability_score(P, W)
) %>% round(4)

## ----session, include=TRUE----------------------------------------------------
sessionInfo()

