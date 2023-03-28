// STAN model for time-to-event outcome
//      assuming exponential baseline hazard
//
// Kamaryn Tanner  17-Mar-2023
//
//  Please see https://ermeel86.github.io/case_studies/surv_stan_example.html
//   for a worked example of a survival model in Stan 



//
// Data Block
//
//    N is the number of observations
//    NC is the number of covariates
//    X is the N by NC matrix of covariates
//    There are N_uncensored + N_censored = N observations
//    times is the N-vector of survival times which is split
//        into 2 vectors, times_censored and times_uncensored
//    There are 2 NC-vectors containing the prior means and standard
//        deviations of the betas (log hazard ratios)
//    We also specify a prior mean and sd for the log baseline hazard (log_lambda)
//
data {
  int<lower=1> N;
  int<lower=1> N_uncensored;
  int<lower=1> N_censored;
  int<lower=0> NC;
  matrix[N_censored, NC] X_censored;
  matrix[N_uncensored, NC] X_uncensored;
  matrix[N,NC] X;
  vector<lower=0>[N_censored] times_censored;
  vector<lower=0>[N_uncensored] times_uncensored;
  vector[NC] prior_mean_betas;
  vector[NC] prior_sd_betas;
  real prior_mean_log_lambda;
  real prior_sd_log_lambda;
}

//
// Parameters block
//
//   parameters are the log hazard ratios (betas) and the 
//       baseline hazard (log_lambda)
parameters {
  vector[NC] betas;
  real log_lambda;
}



//
// Model block
//
//   priors are normally distributed
//   for efficiency, we split the censored and uncensored observations
//   more information on Stan's functions for the exponential distribution
//        can be found at: https://mc-stan.org/docs/functions-reference/exponential-distribution.html
//

model {
  betas ~ normal(prior_mean_betas,prior_sd_betas);
  log_lambda ~ normal(prior_mean_log_lambda,prior_sd_log_lambda);
  target += exponential_lpdf(times_uncensored | exp(log_lambda+X_uncensored*betas));
  target += exponential_lccdf(times_censored | exp(log_lambda+X_censored*betas));
}




//
// Generated Quantities block
//
generated quantities {

}









