######################################
## Sample code to demonstrate model updating process
##    for refitting, intercept recalibration and
##    Bayesian dynamic upating
##
## The code is demonstrated using the pbc dataset. We assume
##    the first 200 records form the original development dataset
##    and the remaining 218 records are new data used to
##    update the original model.
##
## This code accompanies the manuscript entitled:
##   "Dynamic Updating of Clinical Survival Prediction Models 
##    in a Changing Environment"
##
## 28 March 2023, Kamaryn Tanner, London School of Hygiene and Tropical Medicine
##
#######################################



## Load libraries
## ----

library(survival)
library(rstan)

## Settings used in models
## ----

## For Stan:
##    For execution on a local, multicore CPU with excess RAM
options(mc.cores = parallel::detectCores())
##    To avoid recompilation of unchanged Stan programs
rstan_options(auto_write = TRUE)
##    The stan model is in a separate .stan file which must be loaded
sm <- stan_model(file="expSurvivalModel.stan")
CORES <- 2
CHAINS <- 2  ## low chains and iterations for speed in this example
ITER <- 3000 
BURNIN <- 1000
lambdaF <- 0.8  ##the forgetting factor
set.seed(9628)



## Use the publicly available pbc dataset to demonstrate code
## Note: this is not intended to be a proper analysis of the 
##       pbc dataset
## ----
?pbc
dat <- pbc
## Use the combined outcome transplant or death
dat$status2 <- ifelse(dat$status>0, 1, 0)
## Transform bilirubin measurement
dat$log_bili <- log(dat$bili)
## Simplify by using only 3 covariates
dat2 <- dat[,c("id", "time", "status2", "age", "log_bili", "edema")]


## Use first 200 records for original model; 
datOrig <- dat2[1:200,]
## pretend the rest is "new data" for updating the original model
datUpdate <- dat2[201:418,]

## Fit a Cox PH model to obtain the "original" model
## ----
modelOrig <- coxph(Surv(time, status2)~age+edema+log_bili, 
              data=datOrig)
modelOrig


## Update the original model using the next 218 records

## Method 1: Refit
## ----
##   fit a Cox PH model to the new data; no use of information from original model
updateRefit <- coxph(Surv(time, status2)~age+edema+log_bili, 
              data=datUpdate)
updateRefit


## Method 2: Recalibrate
## ----
##   first, get the linear predictor for the new data based on the original model
datUpdate$lp_modelOrig <- predict(modelOrig, newdata=datUpdate, type="lp")
##   now, fit a Cox PH model with linear predictor as an offset
updateRecal <- coxph(Surv(time, status2)~ offset(lp_modelOrig), 
                        data=datUpdate )
updateRecal
##   no change to coefficients; only the baseline hazard changed


## Method 3: Bayesian dynamic update
## ----
## first, set up data to be passed to stan model
X <- datUpdate[,4:6]  #covariate matrix
time <- datUpdate$time
event <- datUpdate$status2
N_censored <- length(which(datUpdate$status2==0))
N <- nrow(X)

## priors for betas derived from modelOriginal
priors <- list()
priors$location <- coefficients(modelOrig)
n <- length(priors$location)
priors$scale <- vector(length=n)
for(i in 1:n){
  priors$scale[i] <- modelOrig$var[i,i]/lambdaF
}
priors$scale <- sqrt(priors$scale)
##can set an uninformative prior on the rate parameter
prior_int <- list(location=0, scale=5)

##create the stan data
stan_data <- list(N=N, N_uncensored=N-N_censored,
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[event==0,]),
                  X_uncensored=as.matrix(X[event!=0,]), X=X,
                  times_censored=time[event==0],
                  times_uncensored=time[event!=0],
                  NC=ncol(X),
                  prior_mean_log_lambda = prior_int$location,
                  prior_sd_log_lambda = prior_int$scale,
                  prior_mean_betas = priors$location,
                  prior_sd_betas = priors$scale)

updateBayes <- sampling(sm, data=stan_data, chains=CHAINS, iter=ITER, 
                   cores=CORES, warmup=BURNIN)


## Updated models are contained in:
##  updateRefit, updateRecal and updateBayes

