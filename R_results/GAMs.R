install.packages("Rcpp")
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_vignettes = TRUE)

library(mgcv)
library(dplyr)
library(ggplot2)

my_data <- read_csv("02-PROKs_selected_columns_forlmms.csv")
View(my_data)
my_data$size_code <- as.factor(my_data$size_code)
my_data$size_code <- relevel(my_data$size_code, ref = "W")


#check normality for link function
hist(my_data$chao1, main = "Histogram of Chao1 Richness", xlab = "Chao1 Richness")

#compare diff models with diff link functions
gam_identity <- gam(chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4) + offset(log(Total)),
                    method = "REML", family = gaussian(link = "identity"), data = my_data)
gam_log <- gam(  chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4) + offset(log(Total)),
                 method = "REML", family = gaussian(link = "log"), data = my_data)
gam_gamma <- gam(  chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4) + offset(log(Total)),
                   method = "REML", family = Gamma(link = "log"), data = my_data)
gam_model_gamma <- gam(chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4) + offset(log(Total)),
                       family = Gamma(link = "log"), data = my_data, method = "REML")


AIC(gam_identity, gam_log, gam_gamma, gam_model_gamma)
#log has the lowest AIC, so richness over time changes according to a log function?

my_data$size_code <- relevel(my_data$size_code, ref = "W")
my_data

gam_model<- gam(
  chao1 ~ size_code + s(weekn, by = size_code, k = 5) + s(depth, k = 4) + offset(log(Total)),
  family = gaussian(link = "log"),
  data = my_data,
  method = "REML"
)

aggregate(chao1 ~ size_code, data = my_data, FUN = mean)


k.check(gam_model)

summary(gam_model)

plot(gam_model, pages = 1)
gam.check(gam_model)
# the model is more uncertain for larger chao1 values (more spread near the end of x axis)


set.seed(123)
train_idx <- sample(1:nrow(my_data), size = 0.8 * nrow(my_data))
train_data <- my_data[train_idx, ]
test_data <- my_data[-train_idx, ]

gam_model <- gam(
  chao1 ~ size_code + 
    s(weekn, by = size_code, k = 5) + 
    s(depth, k = 4) + 
    offset(log(Total)),
  family = Gamma(link = "log"),
  data = train_data,
  method = "REML"
)

predicted <- predict(gam_model, newdata = test_data, type = "response")
actual <- test_data$chao1
plot(actual, predicted, main = "Observed vs Predicted")
abline(0, 1, col = "red")
cor(actual, predicted)  # Correlation between observed and predicted

summary(gam_model)


gam_model_no_offset <- gam(
  chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4),
  family = gaussian(link = "log"),
  data = my_data,
  method = "REML"
)

summary(gam_model_no_offset)


gam_simple <- gam(
  chao1 ~ size_code + offset(log(Total)),
  family = gaussian(link = "log"),
  data = my_data,
  method = "REML"
)
summary(gam_simple)

gam_model <- gam(
  chao1 ~ size_code + s(weekn, by = size_code, k = 8) + s(depth, k = 4) + offset(log(Total)),
  family = gaussian(link = "log"),
  data = my_data,
  method = "REML"
)

summary(gam_model)



gam_model_with_library <- gam(
  chao1 ~ size_code + Total + s(weekn, by = size_code, k = 8) + s(depth, k = 4),
  family = gaussian(link = "log"),
  data = my_data,
  method = "REML"
)
summary(gam_model_with_library)

new_data <- data.frame(
  size_code = c("W", "L", "S", "SL"),
  weekn = median(my_data$weekn),
  depth = median(my_data$depth),
  Total = mean(my_data$Total)
)
new_data$predicted_chao1 <- exp(predict(gam_model_with_library, new_data, type = "link"))

# Compare to observed means
aggregate(chao1 ~ size_code, data = my_data, FUN = mean)
new_data

plot(residuals(gam_model_with_library) ~ my_data$size_code, main = "Residuals by Size Code")

plot(fitted(gam_model), my_data$chao1, main = "Observed vs. Predicted", xlab = "Fitted Values", ylab = "Observed Values")
abline(0, 1, col = "red")  # Perfect fit line

plot(gam_model_with_library, pages = 1)



######## PLOTTING################
#### FIRST PLOT THE OBSERVED RICHNESS WITH LIBRARY SIZE #######
new_data <- data.frame(
  size_code = "W",  # Replace with desired size_code
  Total = seq(min(my_data$Total), max(my_data$Total), length.out = 100),
  weekn = median(my_data$weekn),
  depth = median(my_data$depth)
)
new_data$log_Total <- log(new_data$Total)
new_data$predicted_chao1 <- exp(predict(gam_model_with_library, newdata = new_data, type = "link"))

# Plot observed data and fitted line
ggplot(my_data, aes(x = log(Total), y = chao1)) +
  geom_point(aes(color = size_code), alpha = 0.6) +  # Observed data
  geom_line(data = new_data, aes(x = log_Total, y = predicted_chao1), color = "blue", size = 1) +  # Fitted line
  labs(title = "Observed Data with Fitted GAM Line",
       x = "Log(Total)", y = "Richness (Chao1)") +
  theme_minimal()


######### THEN PLOT THE OBSERVED AND PREDICTED RICHNESS ########
# Prepare prediction data
new_data <- data.frame(
  size_code = levels(my_data$size_code),
  Total = mean(my_data$Total),  # Average library size
  weekn = median(my_data$weekn),
  depth = median(my_data$depth)
)
new_data$log_Total <- log(new_data$Total)
new_data$predicted_chao1 <- exp(predict(gam_model_with_library, newdata = new_data, type = "link"))

observed_means <- aggregate(chao1 ~ size_code, data = my_data, FUN = mean)
combined_data <- merge(observed_means, new_data, by = "size_code")

# Plot
ggplot(combined_data, aes(x = size_code)) +
  geom_bar(aes(y = predicted_chao1, fill = "Predicted"), stat = "identity", alpha = 0.7) +
  geom_point(aes(y = chao1, color = "Observed"), size = 3) +
  labs(#title = "Observed and predicted richness by size fraction",
       x = "Size Fraction", y = "Richness (Chao1)", fill = "", color = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Predicted" = "blue")) +
  scale_color_manual(values = c("Observed" = "red"))
