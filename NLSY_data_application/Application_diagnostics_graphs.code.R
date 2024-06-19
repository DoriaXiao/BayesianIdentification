setwd("C:/Users/xingy/OneDrive/2023 Fall/BayesIdentify/GitHub_code")

# Load datasets
GMM_ML_fit_1c_ageApp <- readRDS("NLSY_data_application/Application_results/NLSY_ageApp/GMM_ML_fit_1c_D10N50.rds")
GMM_ML_fit_2c_ageApp <- readRDS("NLSY_data_application/Application_results/NLSY_ageApp/GMM_ML_fit_2c_D10N50.rds")
GMM_ML_fit_3c_ageApp <- readRDS("NLSY_data_application/Application_results/NLSY_ageApp/GMM_ML_fit_3c_D10N50.rds")
GMM_ML_fit_4c_ageApp <- readRDS("NLSY_data_application/Application_results/NLSY_ageApp/GMM_ML_fit_4c_D10N50.rds")

GMM_ML_fit_1c_occasionApp <- readRDS("NLSY_data_application/Application_results/NLSY_occasionApp/GMM_ML_fit_1c_D10N50.rds")
GMM_ML_fit_2c_occasionApp <- readRDS("NLSY_data_application/Application_results/NLSY_occasionApp/GMM_ML_fit_2c_D10N50.rds")
GMM_ML_fit_3c_occasionApp <- readRDS("NLSY_data_application/Application_results/NLSY_occasionApp/GMM_ML_fit_3c_D10N50.rds")
GMM_ML_fit_4c_occasionApp <- readRDS("NLSY_data_application/Application_results/NLSY_occasionApp/GMM_ML_fit_4c_D10N50.rds")


# Load permuted fit results
fit_permuted_3c_ageApp <- readRDS(file = "NLSY_data_application/Application_results/NLSY_ageApp/fit_permuted_3c_D10N50_d.rds")
fit_permuted_3c_occasionApp <- readRDS(file = "NLSY_data_application/Application_results/NLSY_occasionApp/fit_permuted_3c_D10N50_d.rds")

# Rhat diagnostic
source("C:/Users/xingy/OneDrive/2023 Fall/BayesIdentify/GitHub_code/Diagnostics/Diagnostics.source.R")
priors <- list()
priors[[1]] <- fit_permuted_3c_ageApp

# Check Effective sample size (ESS)
priors[[1]][[3]][,"n_eff"][1:21]

# Traceplot of three lambdas by the specified chains
traceData_ESS(priors, 1, ESS_var = "lambda_1", ESS_chain=4)$coda.ESS.chain

traceplot(Data_reordered=traceData(priors, 1, 5000)$data, 
          num_chains=c(1:5), 
          iterations_per_chain=1000)$traceplot.by.chain

# Rhat diagnostics
rhat_results <- lapply(seq_along(priors), 
                       function(i) find_rhat_positions_cross(priors[[i]], small_threshold = 1.05,  i))

traceData(priors, 1, 5000)$traceplot

# Stuck diagnostic

stuck <- lapply(seq_along(priors), 
                function(i) find_stuck_detail(priors[[i]], i, total_iter = 5000, iter_per_chain = 1000, 10))
stuck

# Figure 1
# Load data
CurranLong_nm <- read_csv("NLSY_data_application/CurranLong_nm.csv")
head(CurranLong_nm)

# Subset data for specified age values
x_values <- c(6:14)
subset <- CurranLong_nm[CurranLong_nm$kidagetv %in% x_values, ]

# Load permuted fit results for ageApp
fit_permuted_3c_D10N50 <- fit_permuted_3c_ageApp

# Define x values
x <- seq(1:5000)

# Extract results and display key statistics
priors[[1]] <- fit_permuted_3c_D10N50
results <- priors[[1]][[3]] %>% round(., 3) %>% as.matrix()
results[, c("mean", "sd", "2.5%", "97.5%", "Rhat")]

# Load necessary scripts
Data <- traceData(priors, 1, 5000)$data

# Extract mean values for lambdas
lambda_means <- c(Data$lambda_1 %>% mean(), Data$lambda_2 %>% mean(), Data$lambda_3 %>% mean())
round(lambda_means, 3)

# Generate sequence for x values
xseq <- seq(0, 8, length.out = 1000)

# Generate y values for each class
yseq3 <- Data$mu_intercept_3 %>% mean() + 
  Data$mu_slope_3 %>% mean() * xseq + 
  Data$sq_mu_slope_3 %>% mean() * xseq^2

yseq1 <- Data$mu_intercept_1 %>% mean() + 
  Data$mu_slope_1 %>% mean() * xseq + 
  Data$sq_mu_slope_1 %>% mean() * xseq^2

yseq2 <- Data$mu_intercept_2 %>% mean() + 
  Data$mu_slope_2 %>% mean() * xseq + 
  Data$sq_mu_slope_2 %>% mean() * xseq^2

# Calculate standard deviations for each class
sd_values3 <- sqrt((Data$sigma1_3 %>% mean())^2 + 
                     (Data$sigma2_3 %>% mean())^2 * xseq^2 + 
                     (Data$sigma1_3 %>% mean() * Data$sigma2_3 %>% mean() * Data$Omega12_3 %>% mean()) * xseq)

sd_values2 <- sqrt((Data$sigma1_2 %>% mean())^2 + 
                     (Data$sigma2_2 %>% mean())^2 * xseq^2 + 
                     (Data$sigma1_2 %>% mean() * Data$sigma2_2 %>% mean() * Data$Omega12_2 %>% mean()) * xseq)

sd_values1 <- sqrt((Data$sigma1_1 %>% mean())^2 + 
                     (Data$sigma2_1 %>% mean())^2 * xseq^2 + 
                     (Data$sigma1_1 %>% mean() * Data$sigma2_1 %>% mean() * Data$Omega12_1 %>% mean()) * xseq)

# Create a subset for age values
subset <- CurranLong_nm[CurranLong_nm$kidagetv %in% x_values, ]

# Plot results
ggplot() +
  geom_boxplot(data = subset, aes(x = factor(kidagetv), y = read), position = position_nudge(x = 0)) +
  geom_line(data = data.frame(x = xseq + 1 , y = yseq1),
            aes(x = x, y = y, color = "Class 1"), size = 1) +
  geom_ribbon(data = data.frame(x = xseq + 1 , y = yseq1, sd_upper = yseq1 + 0.6745 * sd_values1, sd_lower = yseq1 - 0.6745 * sd_values1),
              aes(x = x, ymin = sd_lower, ymax = sd_upper),
              fill = "black", alpha = 0.4) +
  geom_line(data = data.frame(x = xseq + 1 , y = yseq2),
            aes(x = x, y = y, color = "Class 2"), size = 1) +
  geom_ribbon(data = data.frame(x = xseq + 1 , y = yseq2, sd_upper = yseq2 + 0.6745 * sd_values2, sd_lower = yseq2 - 0.6745 * sd_values2),
              aes(x = x, ymin = sd_lower, ymax = sd_upper),
              fill = "mediumblue", alpha = 0.5) +
  geom_line(data = data.frame(x = xseq + 1 , y = yseq3),
            aes(x = x, y = y, color = "Class 3"), size = 1) +
  geom_ribbon(data = data.frame(x = xseq + 1 , y = yseq3, sd_upper = yseq3 + 0.6745 * sd_values3, sd_lower = yseq3 - 0.6745 * sd_values3),
              aes(x = x, ymin = sd_lower, ymax = sd_upper),
              fill = "darkred", alpha = 0.5) +
  scale_y_continuous(breaks = seq(2, 8, by = 1)) +
  labs(x = "Child age",
       y = "Reading recognition",
       color = "Mean trajectories") +
  scale_color_manual(values = c("Class 1" = "black", "Class 2" = "mediumblue", "Class 3" = "darkred")) +
  theme_classic() +
  my_custom_theme() 
           