####---- Libraries ----####
library(TESS) 
library(ape)

####---- Tree ----####

trfn <- file.path( "Datasets/Dataset_S10.tre")
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


### --- CONTINUOUSLY VARIABLE RATE --- ###

prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priorsDecrBD <- c("turnover"=prior_delta,
                  "initial speciation"=prior_lambda,
                  "speciation decay"=prior_alpha)

likelihoodDecrBD <- function(params) {
  speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
  extinction <- function(t) params[1]
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

survivalProbability <- 0.2

prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_time <- function(x) { dunif(x,min=max(times)/2,max=max(times),log=TRUE)}
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
  samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
                             priors = priorsDecrBD,
                             parameters = runif(3,0,1),
                             logTransforms = c(TRUE,TRUE,TRUE),
                             delta = c(1,1,1),
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

marginalLikelihoodDecrBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodDecrBD,
  priors = priorsDecrBD,
  parameters = runif(3,0,1),
  logTransforms = c(TRUE,TRUE,TRUE),
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
                     "DecrBD"=marginalLikelihoodDecrBD,
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

####---- Model Selection with CoMET ----####

# Priors
# number of expected mass extinctions
numExpectedMassExtinctions <- 1

# number of expected rate shifts
numExpectedRateChanges <- 2 # based on BAMM

# Specify the mean and standard deviation of the lognormal
# prior on the speciation rate in real space
speciationPriorMu <- 0.2
speciationPriorSigma <- 0.5

# Specify the mean and standard deviation of the lognormal
# prior on the extinction rate in real space
extinctionPriorMu <- 0.15
extinctionPriorSigma <- 0.5

# Transform the priors on the speciation rate into log space.
speciationRatePriorMean <- log((speciationPriorMu^2)/sqrt(speciationPriorSigma^2+ speciationPriorMu^2))
speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2 /(speciationPriorMu^2)))

# Transform the priors on the extinction rate into log space.
extinctionRatePriorMean <- log((extinctionPriorMu^2) /sqrt(extinctionPriorSigma^2+ extinctionPriorMu^2))
extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2 /(extinctionPriorMu^2)))

# Survival probability in Mass extinction
expectedSurvivalProbability <- 0.2
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *expectedSurvivalProbability /(expectedSurvivalProbability - 1)

# Run CoMET
setwd("Outputs")
set.seed(666)
tess.analysis(tr,
              empiricalHyperPriors = FALSE,
              initialSpeciationRate = speciationPriorMu,
              speciationRatePriorMean = speciationRatePriorMean,
              speciationRatePriorStDev = speciationRatePriorStDev,
              initialExtinctionRate = extinctionPriorMu,
              extinctionRatePriorMean = extinctionRatePriorMean,
              extinctionRatePriorStDev = extinctionRatePriorStDev,
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges,
              numExpectedMassExtinctions = numExpectedMassExtinctions,
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
              MAX_ITERATIONS = 10000,
              dir = "tess_analysis")

# Process CoMET output
output <- tess.process.output("tess_analysis",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

####---- Plot CoMET ----####

# layoput matrix
layout.mat <- matrix(1:9,nrow=3,ncol=3,byrow=TRUE)
layout(layout.mat)

# plot
tess.plot.output(output,fig.types = c("speciation Bayes factors",
                               "speciation shift times",
                               "speciation rates",
                               "extinction Bayes factors",
                               "extinction shift times",
                               "extinction rates",
                               "mass extinction Bayes factors",
                               "mass extinction times",
                               "net-diversification rates"),las=2)

