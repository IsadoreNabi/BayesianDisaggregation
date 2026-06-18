// =============================================================================
// Evidence-based Bayesian disaggregation: state-space model (Model A).
//
// Disaggregates an observed aggregate index (CPI) into K latent sectoral price
// indices phi_{t,k}. Unlike the deterministic family retained as a baseline, the
// aggregate index enters here as GENUINE evidence (an observation density), and
// the cross-sectional concentration and the temporal pattern are real, estimable
// quantities (no cancellation on renormalization). This file is the single
// source of truth for the canonical disaggregation engine.
//
//   Transition (temporal dynamics, partial pooling):
//     log phi_{t,k} = log phi_{t-1,k} + delta_k + tau_k * eta_{t,k},
//                     eta_{t,k} ~ Normal(0, 1),
//     with delta_k ~ Normal(delta_mu, delta_sigma) and
//          log tau_k ~ Normal(log_tau_mu, log_tau_sigma)   (partial pooling).
//
//   Cross-section at t = 1 (estimable concentration, "VAB-anchored" level):
//     log phi_{1,k} = log(phi1_center) + omega_struct * z1_k,  z1_k ~ Normal(0,1),
//     with omega_struct > 0 an ESTIMATED cross-sectional dispersion (the real
//     concentration kappa_struct = 1 / omega_struct^2). This is what the old
//     Dirichlet gamma failed to be: a concentration that does not cancel.
//
//   Observation (CPI as genuine evidence):
//     cpi_t ~ Student-t(nu, sum_k W_{t,k} * phi_{t,k}, sigma_cpi),   [student_obs=1]
//     cpi_t ~ Normal(sum_k W_{t,k} * phi_{t,k}, sigma_cpi),         [student_obs=0]
//   where W are the (known) VAB aggregation weights. The aggregate is strongly
//   identified by this density; the sectoral differentiation is identified by the
//   cross-sectional prior (omega_struct) plus temporal smoothness (tau_k), so the
//   per-sector intervals are honestly wide and prior-influenced. See the vignette.
//
// Non-centered parametrization throughout (z_* are standard-normal primitives).
// Generated quantities expose phi (levels) for downstream use and a pointwise
// log-likelihood for PSIS-LOO.
// =============================================================================

data {
  int<lower=2> T;                       // number of periods
  int<lower=1> K;                       // number of sectors
  vector<lower=0>[T] cpi;               // observed aggregate index (level)
  matrix<lower=0>[T, K] W;              // VAB aggregation weights (rows ~ sum to 1)

  real<lower=0> phi1_center;            // anchor for the initial cross-section level
  real<lower=0> s_omega_struct;         // half-normal scale for omega_struct (log scale)
  real<lower=0> s_delta_mu;             // prior sd of the common drift
  real<lower=0> s_delta_sigma;          // half-normal scale of the drift dispersion
  real          log_tau_loc;            // prior location of log_tau_mu
  real<lower=0> s_log_tau_mu;           // prior sd of log_tau_mu
  real<lower=0> s_log_tau_sigma;        // half-normal scale of log_tau dispersion
  real<lower=0> s_sigma_cpi;            // half-normal scale of the observation sd

  int<lower=0, upper=1> student_obs;    // 1 = Student-t observation, 0 = Gaussian
}

parameters {
  // initial cross-section (non-centered) and its estimable dispersion
  vector[K] z_phi1;
  real<lower=0> omega_struct;

  // drift, partial pooling
  real delta_mu;
  real<lower=0> delta_sigma;
  vector[K] z_delta;

  // innovation scale, partial pooling (log-normal hierarchy)
  real log_tau_mu;
  real<lower=0> log_tau_sigma;
  vector[K] z_tau;

  // transition innovations (non-centered)
  matrix[T - 1, K] z_eta;

  // observation scale and (optional) Student-t degrees of freedom
  real<lower=0> sigma_cpi;
  array[student_obs] real<lower=0> nu_tilde;   // nu = 2 + nu_tilde when student_obs
}

transformed parameters {
  vector[K] delta = delta_mu + delta_sigma * z_delta;
  vector<lower=0>[K] tau = exp(log_tau_mu + log_tau_sigma * z_tau);
  matrix[T, K] log_phi;
  vector[T] agg;                         // VAB-weighted aggregate of phi

  for (k in 1:K)
    log_phi[1, k] = log(phi1_center) + omega_struct * z_phi1[k];
  for (t in 2:T)
    for (k in 1:K)
      log_phi[t, k] = log_phi[t - 1, k] + delta[k] + tau[k] * z_eta[t - 1, k];

  // The exponent is clamped from above before exp() to guard against overflow in
  // pathological warmup proposals (an inf aggregate would reject the step and
  // pollute the chain). The bound 50 (phi ~ e^50) is astronomically above any
  // plausible index level (log phi ~ log 100 ~ 4.6), so it never binds in the
  // stable region: same numerical-hardening logic as bayesianOU D-IMPL-2/9.
  for (t in 1:T) {
    vector[K] phi_t;
    for (k in 1:K) phi_t[k] = exp(fmin(log_phi[t, k], 50));
    agg[t] = W[t] * phi_t;
  }
}

model {
  // ---- priors -----------------------------------------------------------------
  z_phi1        ~ std_normal();
  omega_struct  ~ normal(0, s_omega_struct);        // half-normal (lower bound 0)

  delta_mu      ~ normal(0, s_delta_mu);
  delta_sigma   ~ normal(0, s_delta_sigma);         // half-normal
  z_delta       ~ std_normal();

  log_tau_mu    ~ normal(log_tau_loc, s_log_tau_mu);
  log_tau_sigma ~ normal(0, s_log_tau_sigma);       // half-normal
  z_tau         ~ std_normal();

  to_vector(z_eta) ~ std_normal();

  sigma_cpi     ~ normal(0, s_sigma_cpi);           // half-normal
  if (student_obs)
    nu_tilde[1] ~ gamma(2, 0.1);                    // df concentrated away from 0

  // ---- observation (CPI as genuine evidence) ---------------------------------
  // target += (rather than ~) so the two mutually exclusive branches do not look
  // like a double-update of cpi to the static analyser.
  if (student_obs)
    target += student_t_lpdf(cpi | 2 + nu_tilde[1], agg, sigma_cpi);
  else
    target += normal_lpdf(cpi | agg, sigma_cpi);
}

generated quantities {
  matrix[T, K] phi = exp(log_phi);       // sectoral price-index levels
  real nu = student_obs ? 2 + nu_tilde[1] : 0;   // 0 when Gaussian observation
  vector[T] log_lik;
  {
    real nu_use = student_obs ? 2 + nu_tilde[1] : 0;
    for (t in 1:T)
      log_lik[t] = student_obs
        ? student_t_lpdf(cpi[t] | nu_use, agg[t], sigma_cpi)
        : normal_lpdf(cpi[t] | agg[t], sigma_cpi);
  }
}
