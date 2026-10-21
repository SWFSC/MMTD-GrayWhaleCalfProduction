// =====================================================================
// Fixed-year Negative Binomial model with AR(1) weekly effects  (sensitivity check)
// =====================================================================
// Same as the fixed-year model (GWCAlfCounts_nb_fixedyear.stan)
// EXCEPT the shared weekly shape effects follow an AR(1) process instead of 
// being independent, sum-to-zero draws. This relaxes the "each week is 
// independent" assumption by letting adjacent weeks be correlated 
// (correlation rho, shared across all years).
//
//   rho = 0  recovers the independent model (AR(1) is a strict generalization).
//   rho > 0  smooths the seasonal curve: adjacent weekly effects are correlated.
//
// PURPOSE: sensitivity analysis in response to a reviewer's question about 
// independent weeks in the fixed-year model when weekly counts are obviously
// correlated. 
//
// Expectation (same logic as the shifting-peak model): Because the estimand,
// which is the total number of mother-calf pairs, is the area under the seasonal
// curve, smoothing the weekly deviations redistributes shape without changing
// area, so the abundance estimates should be essentially unchanged. rho is
// estimated from ~12 lag-1 transitions of the single shared 13-week
// sequence.
//
// IMPLEMENTATION NOTE: an AR(1) process is not sum-to-zero on its own, and the
// level/shape additive ridge (add c to every year_eff, subtract c from every
// week_eff) returns if week_eff is left free. We therefore build the AR(1)
// effects on a free vector and then RE-CENTER (subtract the mean) so that
// sum(week_eff) = 0 by construction and year_eff keeps carrying the level,
// exactly as the sum_to_zero_vector did in the fixed-year model. Otherwise,
// MCMC has a difficult time converging
//
// Tomo Eguchi with help from AI (Claude)
// =====================================================================

data {
  int<lower=1> n_obs;
  int<lower=1> n_years;
  int<lower=1> n_weeks;

  array[n_obs] int<lower=0> count_obs;
  vector<lower=0>[n_obs] effort;
  vector[n_obs] log_offset;
  array[n_obs] int<lower=1, upper=n_years> year_idx;
  array[n_obs] int<lower=1, upper=n_weeks> week_idx;
}

parameters {
  real<lower=0, upper=1> p_obs;
  vector[n_years] year_eff;          // FREE, unpooled (no cross-year shrinkage)
  real<lower=0> sigma_week;          // marginal SD of the weekly effects
  real<lower=-1, upper=1> rho;       // AR(1) autocorrelation (shared across years)
  real<lower=0.01> phi;
  vector[n_weeks] week_raw;          // standardized AR(1) innovations (non-centered)
}

transformed parameters {
  vector[n_weeks] week_eff;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  {
    // ---- build the AR(1) sequence (stationary) from standardized innovations ----
    vector[n_weeks] w_ar;
    w_ar[1] = sigma_week * week_raw[1];                 // marginal SD at t=1
    for (t in 2:n_weeks)
      w_ar[t] = rho * w_ar[t - 1]
                + sigma_week * sqrt(1 - square(rho)) * week_raw[t];

    // ---- re-center so the shape sums to zero (level stays in year_eff) ----
    week_eff = w_ar - mean(w_ar);
  }

  for (i in 1:n_obs) {
    log_mu_true[i] = year_eff[year_idx[i]] + week_eff[week_idx[i]];
    log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
  }
}

model {
  year_eff   ~ normal(0, 5);               // wide, INDEPENDENT: no shrinkage
  week_raw   ~ normal(0, 1);               // implies the AR(1) marginal via the recursion above
  rho        ~ normal(0, 0.5);             // weakly informative, favors modest positive correlation
  sigma_week ~ normal(0, 2);               // half-normal (lower bound not set; keep >=0 via prior mass)
  p_obs      ~ beta(20.6967, 2.58418);     // Beta reparam of the Perryman detection estimate (Tinker's suggestion)
  phi        ~ gamma(2, 0.1);

  for (i in 1:n_obs)
    if (effort[i] > 0)
      count_obs[i] ~ neg_binomial_2_log(log_alpha[i], phi);
}

generated quantities {
  vector[n_obs] log_lik;
  vector[n_obs] mu_true;
  vector[n_obs] p_obs_corr;
  array[n_obs] int count_rep;

  for (i in 1:n_obs) {
    mu_true[i] = exp(log_mu_true[i]);
    if (effort[i] > 0) {
      log_lik[i]    = neg_binomial_2_log_lpmf(count_obs[i] | log_alpha[i], phi);
      p_obs_corr[i] = p_obs * exp(log_offset[i]);
      count_rep[i]  = neg_binomial_2_log_rng(log_alpha[i], phi);
    } else {
      log_lik[i]    = 0;
      p_obs_corr[i] = 0;
      count_rep[i]  = 0;
    }
  }
}
