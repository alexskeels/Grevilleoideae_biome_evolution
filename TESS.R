####---- Libraries ----####
library(TESS) 
library(ape)

####---- Tree ----####

trfn <- file.path( "Dataset_S10.tree")
tr <- read.tree(trfn)

####---- Fit TESS Models ----####

times <- as.numeric( branching.times(tr) )
rateChangeTime <- max( times ) / 2
samplingFraction <- 0.73

### --- CONSTANT RATES --- ###
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsConstBD <- c("diversification"=prior_delta,
                   "turnover"=prior_tau)

likelihoodConstBD <- function(params) {
  speciation <- params[1] + params[2]
  extinction <- params[2]
  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = samplingFraction,
                         samplingStrategy = "uniform",
                         log = TRUE)
  return (lnl)
}


### --- EPISODICALLY VARYING RATES --- ###

prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsEpisodicBD <- c("diversification before"=prior_delta_before,
                      "turnover before"=prior_tau_before,
                      "diversification after"=prior_delta_after,
                      "turnover after"=prior_tau_after)

rateChangeTime <- 45 # 45 Ma

likelihoodEpisodicBD <- function(params) {
  speciation <- c(params[1]+params[2],params[3]+params[4])
  extinction <- c(params[2],params[4])
  lnl <- tess.likelihood.rateshift(times,
                                   lambda = speciation,
                                   mu = extinction,
                                   rateChangeTimesLambda = rateChangeTime,
                                   rateChangeTimesMu = rateChangeTime,
                                   samplingProbability = samplingFraction,
                                   samplingStrategy = "uniform",
                                   log = TRUE)
  return (lnl)
}


### --- MASS EXTINCTION EVENT --- ###

survivalProbability <- 0.25

prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_time <- function(x) { dunif(x,min=30,max=70,log=TRUE)}
priorsMassExtinctionBD <- c("diversification"=prior_delta,
                            "turnover"=prior_tau,
                            "mass-extinction time"=prior_time)


likelihoodMassExtinctionBD <- function(params) {
  speciation <- params[1]+params[2]
  extinction <- params[2]
  time <- params[3]
  lnl <- tess.likelihood(times,
                         lambda = speciation,
                         mu = extinction,
                         massExtinctionTimes = time,
                         massExtinctionSurvivalProbabilities =
                           survivalProbability,
                         samplingProbability = samplingFraction,
                         samplingStrategy = "uniform",
                         log = TRUE)
  return (lnl)
}


# Fit models
fit_models <- FALSE

if(fit_models == TRUE){
  set.seed(12345) # remove this line to obtain a random seed
  samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
                              priors = priorsConstBD,
                              parameters = runif(2,0,1),
                              logTransforms = c(TRUE,TRUE),
                              delta = c(1,1),
                              iterations = 10000,
                              burnin = 1000,
                              thinning = 10,
                              adaptive = TRUE,
                              verbose = TRUE)

  set.seed(12345)
  samplesEpisodicBD <- tess.mcmc(likelihoodFunction = likelihoodEpisodicBD,
                                 priors = priorsEpisodicBD,
                                 parameters = runif(4,0,1),
                                 logTransforms = c(TRUE,TRUE,TRUE,TRUE),
                                 delta = c(1,1,1,1),
                                 iterations = 10000,
                                 burnin = 1000,
                                 thinning = 10,
                                 adaptive = TRUE,
                                 verbose = TRUE)
  
  set.seed(12345)
  samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction =
                                         likelihoodMassExtinctionBD,
                                       priors = priorsMassExtinctionBD,
                                       parameters = c(runif(2,0,1),max(times)*3/4),
                                       logTransforms = c(TRUE,TRUE,FALSE),
                                       delta = c(1,1,1),
                                       iterations = 10000,
                                       burnin = 1000,
                                       thinning = 10,
                                       adaptive = TRUE,
                                       verbose = TRUE)
  
}

####---- Evaluate TESS Models ----####

set.seed(12345)
marginalLikelihoodConstBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodConstBD,
  priors = priorsConstBD,
  parameters = runif(2,0,1),
  logTransforms = c(TRUE,TRUE),
  iterations = 1000,
  burnin = 100,
  K = 50)

marginalLikelihoodEpisodicBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodEpisodicBD,
  priors = priorsEpisodicBD,
  parameters = runif(4,0,1),
  logTransforms = c(TRUE,TRUE,TRUE,TRUE),
  iterations = 1000,
  burnin = 100,
  K = 50)

marginalLikelihoodMassExtinctionBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodMassExtinctionBD,
  priors = priorsMassExtinctionBD,
  parameters = c(runif(2,0,1),max(times)*3/4),
  logTransforms = c(TRUE,TRUE,FALSE),
  iterations = 1000,
  burnin = 100,
  K = 50)

# compare fits
candidateModels <- c("ConstBD"=marginalLikelihoodConstBD,
                     "EpisodicBD"=marginalLikelihoodEpisodicBD,
                     "MassExtinctionBD"=marginalLikelihoodMassExtinctionBD)

# Make all possible combinations of the models.
marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels),
                                      M1=names(candidateModels))

# Add a column that is the 2 ln BF for each pair of models.
marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] -
                                    candidateModels[marginalLikelihoodGrid$M1])

# Sort the comparisons by their 2 ln BF in descending order.
marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF,
                                                       decreasing=TRUE),]
marginalLikelihoodGrid
