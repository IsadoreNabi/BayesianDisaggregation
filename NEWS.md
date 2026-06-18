# BayesianDisaggregation 0.2.1

Robustness fix to the CPI reader (debt D-7.2, settled in the OU-hierarchical
Session 10).

* `read_cpi()` now locates the date and CPI columns with `which()` instead of
  `which.max()`. `which.max()` on an all-`FALSE` `str_detect()` (no column
  matched) silently returns index 1, a false positive that pointed `read_cpi` at
  the first column; it also masked ambiguity by taking only the first of several
  matches. `which()` returns every matching index (or `integer(0)`), so the
  existing `length(col) != 1` guard now errors honestly on both no-match and
  multiple-match. Behaviour on a well-formed CPI file is unchanged (verified
  against the bundled `extdata/CPI.xlsx`: 1960-2023, columns `Year`/`CPI`).

# BayesianDisaggregation 0.2.0

Complete redesign to a genuinely evidence-based method. The aggregate index now
enters the estimation as a real observation density; the sectoral indices come
out as a posterior with credible intervals.

## Two honest Bayesian engines

* `disaggregate_statespace()` — **canonical engine**. A Bayesian state-space
  model: a random-walk-with-drift transition in `log phi` (partial pooling on the
  drift and the innovation scale), an estimable cross-sectional concentration,
  and a Student-t (or Gaussian) observation `cpi_t ~ (nu, sum_k W[t,k] phi[t,k],
  sigma)`. Sampled by HMC (cmdstanr or rstan). Returns a `[T, K, draws]` array of
  posterior draws of `phi` — exactly the multiple-imputation input consumed by
  `bayesianOU::fit_ou_nested_mi()` — plus credible bands and diagnostics.
* `disaggregate_conjugate()` — **closed-form Bayesian baseline**. The exact
  linear-Gaussian posterior (Kalman filter + RTS smoother), MCMC-free, with joint
  posterior draws via the Durbin-Koopman simulation smoother. This is the correct
  realization of the package's original "MCMC-free posterior" aspiration: it does
  condition on the aggregate evidence in closed form.

Helpers: `simulate_disagg()` (the model's own DGP, for recovery and examples),
`align_disagg_inputs()` and `disaggregate_from_files()` (read + align CPI and
VAB-weight Excel files), `disagg_default_priors()`, `disagg_stan_code()`.

## Removed: the entire deterministic legacy (F1–F6 audit)

The 0.1.2 "deterministic Bayesian" family never conditioned on the aggregate CPI
(F1): the posterior was derived from the prior weight matrix alone, several
pieces cancelled on renormalization (Dirichlet concentration F2, temporal pattern
F3), the "efficiency" term was a fixed constant (F4), there were no recovery
tests (F5), and `robust_cor` opportunistically picked the larger correlation
(F6). Because that foundational defect cannot be fixed without turning the method
into the new evidence-based engine, the deterministic blend was retained for one
design cycle as a baseline and then removed entirely (it added nothing the two
Bayesian engines do not do, honestly). Deleted: `bayesian_disaggregate()`,
`posterior_weighted/multiplicative/dirichlet/adaptive()`, `compute_L_from_P()`,
`spread_likelihood()`, `coherence_score()`, `numerical_stability_exp()`,
`temporal_stability()`, `stability_composite()`, `interpretability_score()`,
`run_grid_search()`, `save_results()`, and the `robust_cor/kl_divergence/
total_variation/safe_div` utilities.

## Validation and documentation

* Recovery on synthetic data from the model's own DGP (the aggregate is strongly
  identified; the sectoral split is honestly weakly identified — wide bands).
* Bit-for-bit golden regression of the Stan computation path via
  `generate_quantities` on frozen draws (isolating CSV serialization rounding).
* New vignette `evidence-based-disaggregation` documenting the model, the two
  engines, the F1–F6 rationale and the coupling to the nested OU.
* `DESCRIPTION` no longer claims a "novel/original" contribution (anti-overreach).
