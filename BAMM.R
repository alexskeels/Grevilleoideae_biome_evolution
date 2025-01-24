####---- Libraries ----####

library(ape)
library(BAMMtools)
library(phytools)
library(coda)

####---- Prepare Data ----####

# load tree
tree <- read.tree("Datasets/Dataset_S10.tre")

# sampling fraction of different genera in Grevilleoideae
sampling_table <-data.frame(table(sapply(tree$tip.label, FUN=function(x)strsplit(x, "_")[[1]][1])))
colnames(sampling_table) <- c("genus", "n_sampled")

# estimated number of species in each genus
sampling_table$n_estimated <- NA
sampling_table$n_estimated[which(sampling_table$genus == "Alloxylon")] <- 4
sampling_table$n_estimated[which(sampling_table$genus == "Athertonia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Austromuellera")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Banksia")] <- 170
sampling_table$n_estimated[which(sampling_table$genus == "Bleasdalea")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Brabejum")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Buckinghamia")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Cardwellia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Catalepidia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Darlingia")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Embothrium")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Euplassa")] <- 20
sampling_table$n_estimated[which(sampling_table$genus == "Finschia")] <- 3
sampling_table$n_estimated[which(sampling_table$genus == "Floydia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Gevuina")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Grevillea")] <- 360
sampling_table$n_estimated[which(sampling_table$genus == "Hakea")] <- 150
sampling_table$n_estimated[which(sampling_table$genus == "Helicia")] <- 110
sampling_table$n_estimated[which(sampling_table$genus == "Heliciopsis")] <- 14
sampling_table$n_estimated[which(sampling_table$genus == "Hicksbeachia")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Kermadecia")] <- 8
sampling_table$n_estimated[which(sampling_table$genus == "Knightia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Lambertia")] <- 10
sampling_table$n_estimated[which(sampling_table$genus == "Lasjia")] <- 6
sampling_table$n_estimated[which(sampling_table$genus == "Lomatia")] <- 12
sampling_table$n_estimated[which(sampling_table$genus == "Macadamia")] <- 4
sampling_table$n_estimated[which(sampling_table$genus == "Malagasia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Megahertzia")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Musgravea")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Neorites")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Nothorites")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Opisthiolepis")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Oreocallis")] <- 2
sampling_table$n_estimated[which(sampling_table$genus == "Orites")] <- 9
sampling_table$n_estimated[which(sampling_table$genus == "Panopsis")] <- 13
sampling_table$n_estimated[which(sampling_table$genus == "Sleumerodendron")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Sphalmium")] <- 1
sampling_table$n_estimated[which(sampling_table$genus == "Stenocarpus")] <- 22
sampling_table$n_estimated[which(sampling_table$genus == "Strangea")] <- 3
sampling_table$n_estimated[which(sampling_table$genus == "Telopea")] <- 5
sampling_table$n_estimated[which(sampling_table$genus == "Turrillia")] <- 3
sampling_table$n_estimated[which(sampling_table$genus == "Virotia")] <- 6
sampling_table$n_estimated[which(sampling_table$genus == "Xylomelum")] <- 6

# number sampled / number estimated
sampling_table$sampling_fraction <- sampling_table$n_sampled/sampling_table$n_estimated

# vector of sampling fractions
f <- vector("numeric", length(tree$tip.label))

for(i in 1:nrow(sampling_table)){
  f[which(grepl(sampling_table$genus[i], tree$tip.label))] <- sampling_table$sampling_fraction[i]
}

# save file in format for CLaDS model
write.table(f, file="Outputs/CLaDS_sampling_fraction_vector.txt", col.names = F, row.names = F)

# species level sampling fractions table for BAMM
species_sampling_table <- data.frame(speciesName =tree$tip.label, cladeName=NA, samplingFraction=NA)
for(i in 1:nrow(sampling_table)){
  species_sampling_table$samplingFraction[which(grepl(sampling_table$genus[i], tree$tip.label))] <- sampling_table$sampling_fraction[i]
  species_sampling_table$cladeName[which(grepl(sampling_table$genus[i], tree$tip.label))] <- sampling_table$genus[i]
}
write.table(species_sampling_table, file="Outputs/BAMM_sampling_fraction_table.txt", col.names=T, row.names = F)

# estimate priors
setBAMMpriors(tree, outfile="Outputs/BAMM_priors.txt")

# run "Software/divcontrol.txt" in cmd line - e.g., 
  #cd bamm-2.5.0-Windows\bamm-2.5.0-Windows
  #bamm -c Software/BAMM_divcontrol.txt

####---- Analyse Output ----####

# analyse results with bamm tools
# load events and MCMC chain
edata <- getEventData(tree, eventdata = "Outputs/BAMM_event_data.txt", burnin=0.15)
mcmcout <- read.csv("Outputs/BAMM_mcmc_out.txt", header=T)

# look at chains
plot(mcmcout$logLik ~ mcmcout$generation)
plot(mcmcout$N_shifts ~ mcmcout$generation)
plot(mcmcout$eventRate ~ mcmcout$generation)

# see if chains converegd
burnstart <- floor(0.2 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# effective sample size (should be > 200)
effectiveSize(postburn$N_shifts) #yes
effectiveSize(postburn$logLik) # yes
effectiveSize(postburn$eventRate) # yes

# write posterior
write.csv(postburn, file="Outputs/BAMM_posterior.csv")

# now see how many rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
shift_probs <- summary(edata)

postfile <- "Outputs/BAMM_mcmc_out.txt"

# compare prior and posterior - clearly rejects prior
plotPrior(postburn , expectedNumberOfShifts=2)

# plot rates on tree
plot.bammdata(edata, lwd=1, legend=T)

# get credible sets of rate shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)

# plot crdible shift set
plot.credibleshiftset(css)

# get best
par(mfrow=c(1,1))
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 2,logcolor = F , legend=T, spex="netdiv")
addBAMMshifts(best, cex=2.5)
axisPhylo()
          