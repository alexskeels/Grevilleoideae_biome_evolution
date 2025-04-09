# libraries
library(mgcv)
library(ggplot2)
library(cowplot)

# data
climate_space <- as.data.frame(readxl::read_xlsx("Dataset_S3.xlsx", sheet="Climate_Age"))

# remove no climate_richness cells from analysis
climate_space$climate_age                  <- as.numeric(climate_space$climate_age)
climate_space$climate_richness             <- as.numeric(climate_space$climate_richness)
climate_space <- climate_space[-which(climate_space$climate_richness == 0 | is.na(climate_space$climate_richness)),] 

# Fit GAM
gam_model <- gam(log(climate_richness) ~ s(climate_age, k = 4), 
                 data = climate_space, family = gaussian())

# Fit GLM
glm_model <- glm(log(climate_richness) ~ climate_age, 
                    data = climate_space, family = gaussian())

# Compare model fit statistics
AIC_gam <- AIC(gam_model)
AIC_glm <- AIC(glm_model)


# Calculate deviance explained for each model
null_model <- glm(log(climate_richness) ~ 1, 
                  data = climate_space, family = gaussian())
null_dev <- deviance(null_model)

dev_explained_gam <- (null_dev - deviance(gam_model))/null_dev * 100
dev_explained_glm <- (null_dev - deviance(glm_model))/null_dev * 100

# Formal model comparison with ANOVA
anova_result <- anova(glm_model, gam_model, test = "Chisq")

# Create a summary table
comparison <- data.frame(
  Model = c("glm", "GAM with splines"),
  AIC = c(AIC_glm, AIC_gam),
  Deviance_Explained = c(dev_explained_glm, dev_explained_gam)
)

# Print results
print(comparison)
print(anova_result)

# Visual comparison
# Plot the fitted values from both models
plot_data <- data.frame(
  climate_age = climate_space$climate_age,
  climate_richness = climate_space$climate_richness,
  GAM_fitted = predict(gam_model, type = "response"),
  glm_fitted = predict(glm_model, type = "response")
)

# Sort by climate_age for smooth lines
plot_data <- plot_data[order(plot_data$climate_age),]

ggplot(plot_data, aes(x = climate_age, y = log(climate_richness))) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = GAM_fitted, color = "GAM"), size = 1) +
  geom_line(aes(y = glm_fitted, color = "GLM"), size = 1) +
  scale_color_manual(values = c("GAM" = "blue", "GLM" = "red")) +
  labs(title = "",x="climate_age (Ma)",
       y = "log(species climate_richness)", color = "Model") +
  theme_cowplot(14)+
  scale_x_reverse()
