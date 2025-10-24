#' Run grid search for parameter optimization (parallel PSOCK)
#' @param path_cpi Path to CPI Excel file
#' @param path_weights Path to weights Excel file
#' @param grid_df Data frame with parameter combinations (one row = one config)
#' @param n_cores Number of cores for parallel processing
#' @return Data frame with results for all parameter combinations (ordered by composite desc)
#' @export
run_grid_search <- function(path_cpi, path_weights, grid_df,
                            n_cores = parallel::detectCores() - 1) {
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  parallel::clusterEvalQ(cl, {library(BayesianDisaggregation); NULL})
  
  n_cfg <- nrow(grid_df)
  log_msg("INFO", "Running grid search with", n_cfg, "configs on", n_cores, "cores")
  
  results_list <- foreach::foreach(
    i = 1:n_cfg,
    .packages = c("BayesianDisaggregation","readxl","dplyr","tidyr","stringr","purrr"),
    .errorhandling = "pass"
  ) %dopar% {
    cfg <- grid_df[i, ]
    out <- try({
      res <- bayesian_disaggregate(
        path_cpi, path_weights,
        method = cfg$method, lambda = cfg$lambda, gamma = cfg$gamma,
        coh_mult = cfg$coh_mult, coh_const = cfg$coh_const,
        stab_a = cfg$stab_a, stab_b = cfg$stab_b, stab_kappa = cfg$stab_kappa,
        likelihood_pattern = cfg$likelihood_pattern
      )
      as.data.frame(res$metrics)
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) {
      data.frame(
        method = cfg$method, lambda = cfg$lambda, gamma = cfg$gamma,
        coh_mult = cfg$coh_mult, coh_const = cfg$coh_const,
        stab_a = cfg$stab_a, stab_b = cfg$stab_b, stab_kappa = cfg$stab_kappa,
        likelihood_pattern = cfg$likelihood_pattern,
        coherence = NA, stability = NA, interpretability = NA,
        efficiency = NA, composite = NA, T = NA, K = NA
      )
    } else out
  }
  
  dplyr::bind_rows(results_list) |>
    dplyr::mutate(cfg_id = dplyr::row_number()) |>
    dplyr::arrange(dplyr::desc(composite))
}