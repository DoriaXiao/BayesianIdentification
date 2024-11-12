# This script provides an example of log-likelihood calculation outside of Stan.
## Load MCMC samples after post-hoc relabeling
file_paths_rds_D2C5<- list.files("BayesianIdentifation/MCMCResults_D2C5", 
                                    full.names = TRUE,
                                    pattern="rep")  #  Identify RDS files with relabeled MCMC samples
file_contents_rds_D2C5<- file_paths_rds_D2C5%>% 
  map(function(path){
    readRDS(path)
  })

priors <- file_contents_rds_D2C5

## Load Simulated Dataset

data_files<- list.files("BayesianIdentifation/SimDat", 
                        full.names = TRUE,
                        pattern="SimulatedData_a",)  # Identify simulated dataset files
SimDat<- data_files%>% 
  map(function(path){
    read_csv(path)
  })   

# Define time segments and related parameters for calculations
time_seg <- c(0:3)
time_sq_seg <- c(0, 1, 4, 9)
num_seg <- length(time_seg)
num_subjects <- SimDat[[1]]$ID %>% unique() %>% length()  # Number of unique subjects
num_classes <- SimDat[[1]]$latent_group %>% unique() %>% length()  # Number of latent classes

# Organize response data (y) by subject and time segment
y_seg <- array(NA, dim = c(num_subjects, num_seg))
for (j in 1:num_subjects) {
  start_idx <- (j - 1) * num_seg + 1
  end_idx <- j * num_seg
  y_seg[j, ] <- SimDat[[1]]$y[start_idx:end_idx]
}

# Define function for calculating log-likelihood based on mean and covariance
mmn <- function(mu_intercept, mu_slope, sq_mu_slope, time_seg, time_sq_seg, Sigma, sigma_e, y_seg) {
  mu_seg <- mu_intercept + mu_slope * time_seg + sq_mu_slope * time_sq_seg  # Calculate mean vector
  Z_j <- cbind(rep(1, times = length(time_seg)), time_seg)  # Design matrix
  Cov_j <- Z_j %*% Sigma %*% t(Z_j) + diag(sigma_e^2, nrow = length(time_seg))  # Covariance matrix
  
  # Calculate log-likelihood of multivariate normal distribution
  log_likelihood <- dmvnorm(y_seg, mean = mu_seg, sigma = Cov_j, log = TRUE)
  
  return(log_likelihood)
}

# Extract MCMC parameter arrays
file_contents <- file_contents_rds_D2C5[[1]][[2]]
num_iterations <- dim(file_contents)[1]
num_chains <- dim(file_contents)[2]

# Parameter matrices for each MCMC sample
mu_intercept <- file_contents[, , 4:6]
mu_slope <- file_contents[, , 7:9]
sq_mu_slope <- file_contents[, , 10:12]
rho_12 <- file_contents[, , 13:15]
sigma_1 <- file_contents[, , 16:18]
sigma_2 <- file_contents[, , 19:21]
sigma_e <- file_contents[, , 22:24]

# Initialize array to store log-likelihood values
log_likelihood_values <- array(NA, dim = c(num_iterations, num_chains, num_subjects * num_classes))

# Loop through MCMC iterations, chains, subjects, and classes to calculate log-likelihood
for (iter in 1:num_iterations) {
  for (chain in 1:num_chains) {
    for (subject in 1:num_subjects) {
      for (class in 1:num_classes) {
        # Extract parameters for the current class
        mu_intercept_class <- mu_intercept[iter, chain, class]
        mu_slope_class <- mu_slope[iter, chain, class]
        sq_mu_slope_class <- sq_mu_slope[iter, chain, class]
        rho_12_class <- rho_12[iter, chain, class]
        sigma_1_class <- sigma_1[iter, chain, class]
        sigma_2_class <- sigma_2[iter, chain, class]
        sigma_e_class <- sigma_e[iter, chain, class]
        
        # Construct covariance matrix (Sigma) for the current class
        Sigma_class <- matrix(c(sigma_1_class^2,
                                sigma_1_class * sigma_2_class * rho_12_class,
                                sigma_1_class * sigma_2_class * rho_12_class,
                                sigma_2_class^2), nrow = 2, byrow = TRUE)
        
        # Calculate log-likelihood using the mmn function
        log_likelihood_values[iter, chain, (subject - 1) * num_classes + class] <- 
          mmn(mu_intercept_class, mu_slope_class, sq_mu_slope_class, 
              time_seg, time_sq_seg,
              Sigma_class, sigma_e_class, y_seg[subject, ])
      }
    }
  }
}

# Save calculated log-likelihood values to an RDS file
saveRDS(log_likelihood_values, "BayesianIdentifation/MCMCResults_D2C5/log_l.rds")
