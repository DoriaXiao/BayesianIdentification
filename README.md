# Bayesian Identification and Estimation of Growth Mixture Models
This repository provides Stan code for fitting Growth Mixture Models (GMMs) and includes diagnostics and graphs described in our paper (link will be provided upon publication).


## Contents

- Stan Code: Contains Stan code to implement GMMs using fully marginal likelihood and addressing label switching by relabeling techniques.
- Diagnostics for Problematic Behaviors: Includes functions to diagnose stuck sequences, twinlike-class behaviors, and minuscule-class issues.
- Graphs: Provides code to visualize within-class heterogeneity using real data from the NLSY (National Longitudinal Survey of Youth), including class-specific mean trajectories with shaded 50% mid-range and box plots of reading scores. Additionally, it includes graphs visualizing simulation results.

<details>
<summary>Stan Code</summary>

This section includes functions to simulate datasets, compile and run the Stan model, and handle label switching. Detailed code and information are available in the [simulation design code.](Simulation_study/Sim_design.code.R)

### 1. simulate dataset
`data_fun_MCMC` is a function to generate simulated data for MCMC. Refer to the [simulation code](Simulation_study/SimCode.source.R) for more details.
```ruby
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
`Stan_D2C5` is a function to run MCMC with a Dirichlet prior with a concentration parameter of 2 and a Half-Cauchy prior with a scale of 5. Detailed code can be found in the [MCMC code.](Simulation_study/MCMC.source.R)
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

```

### 3. Handling label switching
`pp_sss` is a label switching function that takes a Stan fit object, the number of chains, and the number of iterations as inputs.

</details>
<details>
<summary>Diagnostics for Problematic Behaviors</summary>

### 1. Step 1: Initial Screening based on $\hat{R}$

</details>

