# This file shows an example of log-likelihood calculation outside Stan 

## Load an object after post-hoc relabeling 
file_paths_rds_D2C5<- list.files("BayesianIdentifation/MCMCResults_D2C5", 
                                    full.names = TRUE,
                                    pattern="rep")  # Identify file names
file_contents_rds_D2C5<- file_paths_rds_D2C5%>% 
  map(function(path){
    readRDS(path)
  })

priors <- file_contents_rds_D2C5

## Load Simulated Dataset

data_files<- list.files("BayesianIdentifation/SimDat", 
                        full.names = TRUE,
                        pattern="SimulatedData_a",)  # Identify file names
SimDat<- data_files%>% 
  map(function(path){
    read_csv(path)
  })   

time_seg <- c(0:3)
time_sq_seg <- c(0, 1, 4, 9)
num_seg <- length(time_seg)
num_subjects <- SimDat[[1]]$ID%>%unique()%>%length()
num_classes <- SimDat[[1]]$latent_group%>%unique()%>%length()

y_seg <- array(NA, dim = c(num_subjects, num_seg))

for (j in 1:num_subjects) {
  for (index in 1:num_seg) {
    start_idx <- (j - 1) * num_seg + 1
    end_idx <- j * num_seg
    y_seg[j, ] <- SimDat[[1]]$y[start_idx:end_idx]
  }
}

mmn <- function(mu_intercept, mu_slope, sq_mu_slope, time_seg, time_sq_seg, Sigma, sigma_e, y_seg) {
  mu_seg <- mu_intercept + mu_slope * time_seg + sq_mu_slope * time_sq_seg
  Z_j <- cbind(rep(1, times = length(time_seg)), time_seg)
  Cov_j <- Z_j %*% Sigma %*% t(Z_j) + diag(sigma_e^2, nrow = length(time_seg))
  
  # Calculate the density function using dmnorm
  log_likelihood <- dmvnorm(y_seg, mean = mu_seg, sigma = Cov_j, log = T)
  
  return(log_likelihood)
}

file_contents <- file_contents_rds_D2C5[[1]][[2]]
num_iterations <- dim(file_contents)[1]
num_chains <- dim(file_contents)[2]

mu_intercept <- file_contents[, , 4:6]
mu_slope <- file_contents[, , 7:9]
sq_mu_slope <- file_contents[, , 10:12]
rho_12 <- file_contents[, , 13:15]
sigma_1 <- file_contents[, , 16:18]
sigma_2 <- file_contents[, , 19:21]
sigma_e <- file_contents[, , 22:24]


log_likelihood_values <- array(NA, dim = c(num_iterations, num_chains, num_subjects*num_classes))

# Loop over the specified ranges
for (iter in 1:num_iterations) {
  for (chain in 1:num_chains) {
    for (subject in 1:num_subjects) {
      for (class in 1:num_classes) {
        mu_intercept_class <- mu_intercept[iter, chain, class]
        mu_slope_class <- mu_slope[iter, chain, class]
        sq_mu_slope_class <- sq_mu_slope[iter, chain, class]
        rho_12_class <- rho_12[iter, chain, class]
        sigma_1_class <- sigma_1[iter, chain, class]
        sigma_2_class <- sigma_2[iter, chain, class]
        sigma_e_class <- sigma_e[iter, chain, class]
        
        # Calculate Sigma
        Sigma_class <- matrix(c(sigma_1_class^2,
                                sigma_1_class * sigma_2_class * rho_12_class,
                                sigma_1_class * sigma_2_class * rho_12_class,
                                sigma_2_class^2), nrow = 2, byrow = TRUE)
        
        # Calculate log-likelihood using the mmn function
        log_likelihood_values[iter, chain, (subject - 1) * num_classes + class] <- 
          mmn(mu_intercept_class, mu_slope_class, sq_mu_slope_class, 
              time_seg, time_sq_seg,
              #time_seg[subject, ], time_sq_seg[subject, ],
              Sigma_class, sigma_e_class, y_seg[subject, ])
      }
    }
  }
}

saveRDS(log_likelihood_values, "BayesianIdentifation/MCMCResults_D2C5/log_l.rds")
