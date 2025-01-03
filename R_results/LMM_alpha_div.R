install.packages("Rcpp")
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_vignettes = TRUE)

library(lme4)
library(readr)
library(glmmTMB)
library(ggplot2)
library(emmeans)


your_data <- read_csv("02-PROKs_selected_columns_forlmms.csv")
your_data

your_data$size_code <- as.factor(your_data$size_code)
your_data$depth <- as.factor(your_data$depth)
your_data$weekn <- as.numeric(your_data$weekn)
your_data$nASVs <- as.numeric(your_data$nASVs)
your_data$Total <- as.numeric(your_data$Total)

lm_model <- lm(nASVs ~ size_code * depth + weekn, data = your_data)
summary(lm_model)

test_model <- glmmTMB(nASVs ~ size_code * depth + weekn + (1 | depth), data = your_data)
summary(test_model)

your_data$predicted <- predict(test_model, type = "response")

ggplot(your_data, aes(x = weekn, y = nASVs, color = depth)) +
  geom_point() +  # Observed data
  geom_line(aes(y = predicted)) +  # Model predictions
  facet_wrap(~ size_code) +  # Separate by size code
  theme_minimal() +
  labs(title = "Temporal Trends in Richness",
       x = "Week",
       y = "Richness (nASVs)")


ggplot(your_data, aes(x = weekn, y = predicted, color = depth)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Predicted Richness Trends by Depth",
       x = "Week",
       y = "Predicted Richness (nASVs)")

# Post-hoc pairwise comparisons
emmeans_results <- emmeans(test_model, pairwise ~ size_code | depth)
summary(emmeans_results)
plot(emmeans_results)

#---------------------------------depth-specific---------------------------------------
depth_specific_data <- subset(your_data, depth == 1)

model_depth_specific <- lmer(nASVs ~ size_code + weekn + (1 | weekn),
                             data = depth_specific_data)
summary(model_depth_specific)

# Marginal means for size_code
emmeans_size_code <- emmeans(model_depth_specific, ~ size_code)
summary(emmeans_size_code)

# Pairwise comparisons between size fractions
emmeans_pairwise <- emmeans(model_depth_specific, pairwise ~ size_code)
summary(emmeans_pairwise)




depth_levels <- unique(your_data$depth)
results_list <- list()

for (depth in depth_levels) {
  depth_data <- subset(your_data, depth == depth)
  model <- lmer(nASVs ~ size_code + weekn + (1 | weekn), data = depth_data)
  results_list[[depth]] <- emmeans(model, ~ size_code)
}



results_df <- do.call(rbind, lapply(names(results_list), function(depth) {
  data.frame(Depth = depth, as.data.frame(results_list[[depth]]))
}))

# Plot results
library(ggplot2)

ggplot(results_df, aes(x = size_code, y = emmean, color = Depth)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(title = "Estimated Richness by Size Fraction at Each Depth",
       x = "Size Fraction",
       y = "Estimated Richness (nASVs)") +
  theme_minimal()


##############################################DTW#############################
install.packages("dtw")
install.packages("vegan")
install.packages("proxy")

library(vegan)
library(reshape2)

# Subset data for a specific depth
df_depth <- subset(your_data, depth == 10) # Example for depth 10

# Reshape data into size_code x weekn matrix
data_matrix <- dcast(df_depth, size_code ~ weekn, value.var = "nASVs")
rownames(data_matrix) <- data_matrix$size_code
data_matrix <- data_matrix[-1]  # Remove size_code column
data_matrix

library(dtw)

# Pairwise DTW distances
dtw_dist <- dist(data_matrix, method = function(x, y) dtw(x, y, step.pattern = symmetricP2)$distance)

# Hierarchical clustering based on DTW distances
hclust_res <- hclust(dtw_dist, method = "average")
plot(hclust_res, main = "DTW Dendrogram of Size Fractions (Depth 10)")







##########################################################MODEL#########
library(mgcv)
model <- gam(nASVs ~ s(weekn, bs = "cs") + Total + size_code + depth,
             family = poisson(link = "log"),  # Adjust to "negative binomial" if overdispersion
             data = your_data)
gam.check(model)
summary(model)

# Create a prediction data frame
new_data <- expand.grid(weekn = seq(min(your_data$weekn), max(your_data$weekn), length.out = 100),
                        Total = mean(your_data$Total),
                        size_code = unique(your_data$size_code),
                        depth = unique(your_data$depth))

new_data$PredictedRichness <- predict(model, new_data, type = "response")

# Plot
ggplot(new_data, aes(x = weekn, y = PredictedRichness, color = size_code)) +
  geom_line() +
  facet_wrap(~ depth) +
  labs(title = "Temporal Trends in Richness by Size Fraction",
       x = "Time (Weeks)", y = "Predicted Richness") +
  theme_minimal()




overdispersion <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
overdispersion

library(mgcv)
your_data$size_code <- as.factor(your_data$size_code)
your_data$size_code <- relevel(your_data$size_code, ref = "W")
model_nb <- gam(nASVs ~ s(weekn, bs = "cs") + Total + size_code + depth,
                family = nb(),  # Negative binomial distribution
                data = your_data)
summary(model_nb)

gam.check(model_nb)

AIC(model, model_nb)

# Generate predictions
new_data <- data.frame(weekn = seq(min(your_data$weekn), max(your_data$weekn), length.out = 100),
                       Total = mean(your_data$Total),
                       size_code = "L",
                       depth = mean(your_data$depth))

new_data$nASVs <- predict(model_nb, newdata = new_data, type = "response")

# Plot
ggplot(new_data, aes(x = weekn, y = nASVs)) +
  geom_line() +
  labs(title = "Temporal Trend in Richness (Negative Binomial GAM)", 
       x = "Week", y = "Predicted Richness") +
  theme_minimal()


# Fit GLMMTMB model with negative binomial family
glmm_nb <- glmmTMB(nASVs ~ weekn + Total + size_code + depth + (1 | depth),
                   family = nbinom2,  # Negative binomial with quadratic variance
                   data = your_data)
summary(glmm_nb)

AIC(model, model_nb, glmm_nb)


# Generate prediction data
new_data <- expand.grid(weekn = seq(min(your_data$weekn), max(your_data$weekn), length.out = 100),
                        Total = mean(your_data$Total),  # Average sequencing depth
                        size_code = unique(your_data$size_code),  # All size fractions
                        depth = unique(your_data$depth))  # All depths

# Add predicted richness
new_data$nASVs <- predict(model_nb, newdata = new_data, type = "response")

ggplot(new_data, aes(x = weekn, y = nASVs, color = size_code)) +
  geom_line(size = 1) +
  geom_point(data = your_data, aes(x = weekn, y = nASVs, color = size_code), alpha = 0.5) +
  facet_wrap(~ depth, scales = "free_y") +
  labs(title = "Temporal Trends in Richness by Depth and Size Fraction",
       x = "Week (Time)", 
       y = "Predicted Richness",
       color = "Size Fraction") +
  theme_minimal()

summary(model_nb)


new_data <- data.frame(size_code = unique(your_data$size_code),
                       Total = mean(your_data$Total),
                       depth = mean(your_data$depth),
                       weekn = mean(your_data$weekn))
new_data$nASVs <- predict(model_nb, newdata = new_data, type = "response")

ggplot(new_data, aes(x = size_code, y = nASVs, fill = size_code)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Predicted Richness by Size Fraction",
       x = "Size Fraction", y = "Predicted Richness") +
  theme_minimal()

