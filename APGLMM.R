####---- Libraries ----####

library(ape)
library(phytools)
library(BioGeoBEARS)
library(tidyverse)
library(ggtree)
library(treeio)
library(ltstR)
library(plyr)
library(reshape2)

# IMPORTANT: Load Dr. Russell Dinnage's branch of phyr which conatins the ancestral PGLMM
if(require("phyr", quietly = TRUE)) {
  ver <- packageVersion("phyr")
  if(ver != "1.1.9") {
    devtools::install_github("daijiang/phyr", ref = "ancestral")
  }
} else {
  devtools::install_github("daijiang/phyr", ref = "ancestral")
}

library(MCMCglmm)
library(INLA)
library(ggplot2)
library(phyr)

####---- Functions ----####

# Bayes factors
BF <- function(mod_1, mod_2, prior_p = c(0.5, 0.5)) {
  bf <- (mod_1$logLik - mod_2$logLik) * (prior_p[1] / prior_p[2])
  bf
}

# get most likely state of each node in the tree
stateProbs <- function(results_object, trfn, geogfn, statenames, maxrange=3){
  # get tipranges
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  tipranges
  
  areas <- getareas_from_tipranges_object(tipranges)
  
  numstates = cladoRcpp::numstates_from_numareas(numareas = length(areas), 
                                                 maxareas = maxrange, include_null_range = results_object$inputs$include_null_range)
  
  states_list_areaLetters = areas_list_to_states_list_new(areas, 
                                                          maxareas =maxrange, include_null_range = results_object$inputs$include_null_range)
  states_list_0based = cladoRcpp::rcpp_areas_list_to_states_list(areas, 
                                                                 maxareas = maxrange, include_null_range = results_object$inputs$include_null_range)
  
  
  # Make the list of ranges
  ranges_list = NULL
  for (i in 1:length(states_list_0based))
  {    
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
    {
      tmprange = "_"
    } else {
      tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
    }
    ranges_list = c(ranges_list, tmprange)
  }
  
  # Look at the ranges list
  ranges_list
  
  # get probabilities of each state
  probs <- as.data.frame( results_object$ML_marginal_prob_each_state_at_branch_top_AT_node) 
  colnames(probs) <- ranges_list
  
  return(probs)
}

####---- Load Data ----####

# load BGB
load("Datasets/Dataset_S17_DEC.Rdata")

# load CLaDS results from Julia 
load("Outputs/ClaDS_output.RData")

####---- BGB Info ----####

# model
model <- "timestrat_DEC_m_w_n"

# max range size
max_range_size = 3

# get geography
geogfn   <- file.path("Outputs/koppen_regional_geog.txt")
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

# get phylogeny 
tr <- read.tree("Datasets/Dataset_S12_species_level_phylo.tre")

# state names
statenames <- c("-","a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k") # vector of georegion codes
names(statenames) <- 0:11

# probabilities of each state
probs <- stateProbs(res,  trfn, geogfn, statenames)
states <- colnames(probs)[apply(probs, 1, which.max)]

# get diversification rate information by branch
rates <- CladsOutput$lambdai_map

# phylo data frame
tr <- force.ultrametric(tr, method="extend")
tree_data <- prt(tr)

# add states
tree_data$states <- states

####---- Diversification Rates from ClaDS ----####

# remove root (no rate estimated for root)
tree_data <- tree_data[-c(694),] 

# put rates in data frame
tree_data$rates <- NA
for(i in 1:nrow(tree_data)){
  tree_data$rates[i] <- rates[which(tr$edge[,1] == tree_data$ancestor[i] & tr$edge[,2] == tree_data$node[i])]
}
colnames(tree_data)[which(colnames(tree_data)=="label")] <- "BGB_label"

# merge with ggtree / treeio format 
tree_data_tib <- as_tibble(tr)
tree_data <- merge(tree_data_tib, tree_data, by.x=c("parent", "node"), by.y=c("ancestor", "node"), all=T)

# order tree
tree_data <- tree_data[order(tree_data $node, decreasing=F),]
rates_df <- tree_data[, c( "parent",  "node", "branch.length", "label", "states", "rates", "node_ht")]

####---- Binary States ----####

# add a series of binary predictors for the states
rates_df$AustraloPapuan_Tropical <- ifelse(grepl("a", rates_df$states), 1, 0)         
rates_df$AustraloPapuan_Subtropical <- ifelse(grepl("b", rates_df$states), 1, 0)         
rates_df$AustraloPapuan_Mediterranean <- ifelse(grepl("c", rates_df$states), 1, 0)         
rates_df$AustraloPapuan_SemiArid <- ifelse(grepl("d", rates_df$states), 1, 0)         
rates_df$AustraloPapuan_Arid <- ifelse(grepl("e", rates_df$states), 1, 0)         
rates_df$SthAmerica <- ifelse(grepl("f", rates_df$states), 1, 0)         
rates_df$Madagascar <- ifelse(grepl("g", rates_df$states), 1, 0)         
rates_df$Cape <- ifelse(grepl("h", rates_df$states), 1, 0)         
rates_df$TropicalAsia <- ifelse(grepl("i", rates_df$states), 1, 0)         
rates_df$NewCaledonia <- ifelse(grepl("j", rates_df$states), 1, 0)         
rates_df$NewZealand <- ifelse(grepl("k", rates_df$states), 1, 0)      

####---- Time Bins ----####

# Get the time bins to match the resolution of the paleobiome data
rates_df$time <- round(max(rates_df$node_ht, na.rm=T) - rates_df$node_ht, 3)
rates_df$time_bin <- round_any(rates_df$time, 10, ceiling)

####---- Area ----####

# get the area of the biome at each time step
rates_df$ area <- 0
area <- read.csv("Outputs/biome_paleoarea.csv")

area$value[which(is.na(area$value))] <- 0
for(region in unique(colnames(rates_df)[8:18])){
  for(time in unique(rates_df$time_bin)){
    if( region %in% area$variable){
      rates_df$area[which( rates_df[,region]==1 & rates_df$time_bin == time)] <- 
        rates_df$area[which( rates_df[,region]==1 & rates_df$time_bin == time)] + area$value[which(area$variable==region & area$time==time)]}
  }
}

####---- Time since Biome Occupation ----####

rates_df$biome_age <- 0
for(region in unique(colnames(rates_df)[8:18])){
  
  # first colonisation of biome
  biome_origination_date <- max(rates_df$time[which(rates_df[, region]==1)])
  
  
  for(i in 1:length(rates_df$biome_age[which(rates_df[, region]==1)])){
    
    # time since biome was first colonised is larger than time since other occupied biomes were colonised
    if(biome_origination_date - rates_df$time[which(rates_df[, region]==1)][i] > rates_df$biome_age[which(rates_df[, region]==1)][i]){
      rates_df$biome_age[which(rates_df[, region]==1)][i] <- biome_origination_date - rates_df$time[which(rates_df[, region]==1)][i]
    }
  }
}

####---- Biome Shifting ----####
rates_df$biome_shift <- 0

for(i in 1:nrow(rates_df)){
  rates_df$biome_shift[i] <- ifelse(any(!strsplit(rates_df$states[i], "")[[1]] %in% strsplit(rates_df$states[which(rates_df$node == rates_df$parent[i])], "")[[1]]), 1, 0)
}

####---- Standing Diversity ----####

source("Software/LTSTR.R")

# get diversity through time
# make sure no two nodes have exactly the same height
CET <- correctNodeHeights2(tree_data)

# constrain the root
CET$node.type[694] <- "root"
CET$time_bp[694] <- max(branching.times(tr))
CET$states[694] <- "a"

# Use LTSTR package functions to extract diversity through time curves in each state
ETT <- getEventTiming2(tr, CET)
LTST_dt <- getLTSTDataTable2(ETT)
LTST_dt2 <- melt(LTST_dt, id.vars = 'time', variable.name = 'region')
rates_df$diversity <- NA
for(i in 1:nrow(rates_df)){
  
  states <- strsplit(rates_df$states[i], "")[[1]]
  time <- rates_df$time[i]
  
  older_times <- LTST_dt$time[which(LTST_dt$time-time > 0)] # has to be before not after the node
  
  closest_time <- older_times[which(older_times == min(older_times))]
   
  diversity_closest_time <- sum(LTST_dt[which(LTST_dt$time==closest_time), states])
  
  rates_df$diversity[i] <- diversity_closest_time
}

####---- Standing Diversity ----####

# confirm tree
tree <- tr
tree$edge.length[tree$edge.length == 0] <- 0.5
tree <- force.ultrametric(tree, method = "extend")

# do the tree and table match?
all(na.omit(rates_df$label) %in% tree$tip.label)
all(tree$tip.label %in% rates_df$label)

#calculate inverse phylogenetic covariance matrix
tree$edge.length <- tree$edge.length / max(cophenetic(tree))
phy_Ainv <- inverseA(tree, scale = FALSE)$Ainv

# add node labels to tree
tree$node.label <- paste0("Node", seq_len(tree$Nnode))

# create a new dataset with matching node labels
data_phyr <- rates_df 
data_phyr$label[is.na(data_phyr$label)] <- tree$node.label[data_phyr$node[is.na(data_phyr$label)] - tree$Nnode - 1]
data_phyr$phy <- data_phyr$label

## remove root node which has no data and is not in the ancestral phylogenetic covariance
## matrix. phyr doesn't let you use any data with labels not in the covariance, to reduce
## the chance of errors
data_phyr <- data_phyr %>%
  filter(node != 694)

# Na areas should be 0
data_phyr$area[which(is.na(data_phyr$area))] <- 0

# log transform area and diversity
data_phyr$area <- log(data_phyr$area + 1)
data_phyr$diversity <- log(data_phyr$diversity + 1)

# scale and center all varaibles
data_phyr[, c(6,19:24)] <- scale(data_phyr[, c(6,19:24)])

# fit a basic linear model without accoutning for non-independence of data
model_formula_m1 <- formula(log(rates+1) ~  area + biome_age + biome_shift + diversity)
lm1 <- lm(model_formula_m1, data=data_phyr)
summary(lm1 )

## we calculate an unbiased estimate of the brownian motion evolutionary rate
## on model residuals, which we base a weakly informative prior on.
sdres_rates <- pic(resid(lm(model_formula_m1, 
                             data = tibble(label = tree$tip.label) %>%
                               left_join(data_phyr))), tree, scaled = FALSE) 
sdres_rates <- mean(sdres_rates^2)
## we multiply by three to be conservative, our prior will be:
## most of the prior probability density is less than 3 times the
## expected rate under brownian motion, this shrinks estimates
## slightly towards zero and assentially removes overly high
## rates from being estimated

sdres_rates <- sdres_rates * 3

# Fit intercept only model
rates_model_m0 <- pglmm(log(rates+1)  ~ 1 +
                          (1 | phy__), 
                        cov_ranef = list(phy = tree),
                        ancestral = "phy",
                        data = data_phyr,
                        prior = "pc.prior",
                        prior_mu = sdres_rates,
                        prior_alpha = 0.01,
                        verbatim_mode = TRUE,
                        bayes = TRUE)

# Fit model with predictors
rates_model_m1 <- pglmm(log(rates+1) ~ 
                            area +
                            biome_age + 
                            biome_shift +
                            diversity + 
                            (1 | phy__), 
                          cov_ranef = list(phy = tree),
                          ancestral = "phy",
                          data = data_phyr,
                          prior = "pc.prior",
                          prior_mu = sdres_rates,
                          prior_alpha = 0.01,
                          verbatim_mode = TRUE,
                          bayes = TRUE)


# confirm its better than inetrcept only
BF(rates_model_m0, rates_model_m1)

# ok look at coefficients
plot_bayes(rates_model_m1)

# look at model summaries
summary(rates_model_m1)

# R-squared
rr2::R2(rates_model_m0)
rr2::R2(rates_model_m1)