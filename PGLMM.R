####---- Libraries ----####

# IMPORTANT: Load Dr. Russell Dinnage's branch of phyr which conatins the ancestral PGLMM
if(require("phyr", quietly = TRUE)) {
  ver <- packageVersion("phyr")
  if(ver != "1.1.9") {
    devtools::install_github("daijiang/phyr", ref = "ancestral")
  }
} else {
  devtools::install_github("daijiang/phyr", ref = "ancestral")
}
library(phyr)
library(MCMCglmm)
library(INLA)
library(ggplot2)
library(cowplot)
library(ape)
library(dplyr)

####---- Functions ----####

# Bayes factors
BF <- function(mod_1, mod_2, prior_p = c(0.5, 0.5)) {
  bf <- (mod_1$logLik - mod_2$logLik) * (prior_p[1] / prior_p[2])
  bf
}

data_phyr <- readxl::read_xlsx("Dataset_S3.xlsx", sheet="Node Data Table") 

# scale and center
data_phyr$area <- log(data_phyr$area + 1)
data_phyr$sympatric_species_richness <- log(data_phyr$sympatric_species_richness)
data_phyr$rates <- log(data_phyr$rates)
data_phyr$node_tip <- ifelse(data_phyr$node_type == "tip", 1, 0)



data_phyr[, c(6:9, 11)] <- scale(data_phyr[, c(6:9, 11)])

p1 <- ggplot(data_phyr, aes(y=rates, x=as.factor(biome_shift), fill=node_type))+
  geom_boxplot(alpha=0.7)+
  theme_cowplot(14)

p2 <- ggplot(data_phyr, aes(y=rates, x=area))+
  geom_point(alpha=0.7, mapping=aes(colour=node_type))+
  geom_smooth(method="lm")+
  theme_cowplot(14)

p3 <- ggplot(data_phyr, aes(y=rates, x=occupation_time))+
  geom_point(alpha=0.7, mapping=aes(colour=node_type))+
  geom_smooth(method="lm")+
  theme_cowplot(14)

p4 <- ggplot(data_phyr, aes(y=rates, x=sympatric_species_richness))+
  geom_point(alpha=0.7, mapping=aes(colour=node_type))+
  geom_smooth(method="lm")+
  theme_cowplot(14)

cowplot::plot_grid(p1, p2, p3, p4)

## we calculate an unbiased estimate of the brownian motion evolutionary rate
## on model residuals, which we base a weakly informative prior on.
model_formula_m1 <- formula(rates ~   area +
                              occupation_time + 
                              biome_shift+
                              sympatric_species_richness + time_bp + node_tip)
lm1 <- lm(model_formula_m1, data=data_phyr)
summary(lm1 )

model_formula_m0 <- formula(rates ~   1)
lm0 <- lm(model_formula_m0, data=data_phyr)
summary(lm0 )

anova(lm0, lm1)

# read tree in
tree <- read.tree("Dataset_S10.tree")

data_phyr$label[is.na(data_phyr$label)] <- tree$node.label[data_phyr$node[is.na(data_phyr$label)] - tree$Nnode - 1]
data_phyr$phy <- data_phyr$label

## remove root node which has no data and is not in the ancestral phylogenetic covariance
## matrix. phyr doesn't let you use any data with labels not in the covariance, to reduce
## the chance of errors
data_phyr <- data_phyr %>%
  filter(node != 694)

pic_df <- tibble(label = tree$tip.label) %>% 
  left_join(data_phyr)
rownames(pic_df) <- pic_df$label



sdres_rates_m0 <- pic(resid(lm(model_formula_m0, data = pic_df)), tree, scaled = FALSE) 
sdres_rates_m0 <- mean(sdres_rates_m0^2)

sdres_rates_m1 <- pic(resid(lm(model_formula_m1, data = pic_df)), tree, scaled = FALSE) 
sdres_rates_m1 <- mean(sdres_rates_m1^2)

## we multiply by three to be conservative, our prior will be:
## most of the prior probability density is less than 3 times the
## expected rate under brownian motion, this shrinks estimates
## slightly towards zero and assentially removes overly high
## rates from being estimated
sdres_rates_m0 <- sdres_rates_m0 * 3

sdres_rates_m1 <- sdres_rates_m1 * 3

rates_model_m0_1 <- pglmm(rates  ~ 1 +
                            (1 | phy__), 
                          cov_ranef = list(phy = tree),
                          family="gaussian",
                          ancestral = "phy",
                          data = data_phyr,
                          prior = "pc.prior.auto",
                          prior_mu = sdres_rates_m0,
                          prior_alpha = 0.01,
                          verbatim_mode = TRUE,
                          bayes = TRUE)

rates_model_m1_1 <- pglmm(rates ~ 
                            area +
                            occupation_time + 
                            biome_shift +
                            sympatric_species_richness + time_bp + node_tip + (1 | phy__), 
                          cov_ranef = list(phy = tree),
                          ancestral = "phy",
                          data = data_phyr,
                          prior = "pc.prior",
                          prior_mu = sdres_rates_m1,
                          prior_alpha = 0.01,
                          verbatim_mode = TRUE,
                          bayes = TRUE)



# which model is best?
BF(rates_model_m1_1, rates_model_m0_1)

plot_bayes(rates_model_m1_1)


# sensistivity test: presented in Appendix S1

# compare two similar model. 1 with node+tip, once with tip-only. Keep only shared rpedictors.

rates_model_m2 <- pglmm(rates ~ 
                          area +
                          occupation_time + 
                          biome_shift +
                          sympatric_species_richness +  (1 | phy__), 
                        cov_ranef = list(phy = tree),
                        ancestral = "phy",
                        data = data_phyr,
                        prior = "pc.prior",
                        prior_mu = sdres_rates_m1,
                        prior_alpha = 0.01,
                        verbatim_mode = TRUE,
                        bayes = TRUE)

rates_model_m2_tip <- pglmm(rates ~ 
                              area +
                              occupation_time + 
                              biome_shift +
                              sympatric_species_richness + (1 | phy__), 
                            cov_ranef = list(phy = tree),
                            #ancestral = "phy",
                            data = data_phyr[which(data_phyr$node_type=="tip"),],
                            prior = "pc.prior",
                            prior_mu = sdres_rates_m1,
                            prior_alpha = 0.01,
                            verbatim_mode = TRUE,
                            bayes = TRUE)

# fig S12

plot_bayes(rates_model_m2)

plot_bayes(rates_model_m2_tip)
