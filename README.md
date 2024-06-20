# Bayesian Identification and Estimation of Growth Mixture Models
This repository provides Stan code for fitting Growth Mixture Models (GMMs) and includes diagnostics and graphs described in our paper (link will be provided upon publication).


## Contents

- Stan Code: Contains Stan code to implement GMMs using fully marginal likelihood and addressing label switching by relabeling techniques.
- Diagnostics for Problematic Behaviors: Includes functions to diagnose stuck sequences, twinlike-class behaviors, and minuscule-class issues.
- Graphs: Provides code to visualize within-class heterogeneity using real data from the NLSY (National Longitudinal Survey of Youth), including class-specific mean trajectories with shaded 50% mid-range and box plots of reading scores. Additionally, it includes graphs visualizing simulation results.

<details>

<summary>Stan Code</summary>

This section includes functions to simulate datasets, compile and run the Stan model, and handle label switching. Detailed code and information are available in the [simulation design code](Simulation_study/Sim_design.code.R)

### 1. simulate dataset
```ruby
# Source the R script containing the function to simulate data
source("~/SimCode.source.R")

# Function to generate simulated data for MCMC
data_fun_MCMC <- function() {
  # Generate simulated data using the UncGMM_data function
  UncGMM_dat <- UncGMM_data(
    n_pers = 405,
    n_time = 4,
    beta_int = c(3.417, 2.637, 2.019), 
    beta_slo_time = c(1.297, 2.607, 1.393), 
    beta_slo_time_sq = c(-0.094, -0.472, -0.120),
    sd_i = c(0.967, 0.599, 0.188), 
    sd_s = c(0.263, 0.265, 0.356),  
    cor_is = c(-0.567, 0.048, 0.677), 
    sd_r = 0.469, 
    K = 3,
    lambda_K = c(0.266, 0.269, 0.465)
  )
  return(as.data.frame(UncGMM_dat))  # Return simulated data as data frame
}

# Load saved simulated datasets
data_files <- list.files("~/SimDat", 
                         full.names = TRUE,
                         pattern = "SimulatedData_a")  # Identify file names

SimDat <- data_files %>% 
  map(function(path) {
    read_csv(path)  # Read each dataset into a list
  })   
```
### 2. Running model
```ruby
# Install CmdStan with specified number of cores
install_cmdstan(cores = 4)  # Insert your number of cores
cmdstan_path()  # Check the path where CmdStan is installed
cmdstan_version()  # Check the version of CmdStan installed

# Compiling a GMM with marginal likelihood
sq_GMM_ML_file <- file.path(cmdstan_path(), "GMM_ML.stan")
sq_GMM_ML_mod <- cmdstan_model(sq_GMM_ML_file)

# Saving Stan objects for each dataset
for (i in 1:length(SimDat)) {
  source("~/MCMC.source.R")  # Source the R script containing Stan_D2C5 function
  setwd("~/BayesIdentify/MCMCResults_D2C5")  # Set working directory for saving results
  Stan_D2C5(SimDat[[i]], K = 3)$save_object(sprintf("Stan_a%i.rds", i))  # Save Stan objects
}

# Stan_D2C5 function in MCMC.source.R

```

### 3. Handling label switching
`pp_sss` is a label switching function that takes a Stan fit object, the number of chains, and the number of iterations as inputs.

