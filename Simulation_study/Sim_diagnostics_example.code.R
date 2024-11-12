# Load necessary source file for data simulation
source("~/Simulation_study/SimCode.source.R")

# Function to simulate data using specified parameters for an Unconditional Growth Mixture Model (UncGMM)
data_fun_MCMC <- function() {
  UncGMM_dat <- UncGMM_data(
    n_pers = 405,  # Number of individuals
    n_time = 4,  # Number of time points
    beta_int = c(3.417, 2.637, 2.019),  # Intercept values for each class
    beta_slo_time = c(1.297, 2.607, 1.393),  # Linear slope coefficients for each class
    beta_slo_time_sq = c(-0.094, -0.472, -0.120),  # Quadratic slope coefficients for each class
    sd_i = c(0.967, 0.599, 0.188),  # Intercept standard deviations for each class
    sd_s = c(0.263, 0.265, 0.356),  # Slope standard deviations for each class
    cor_is = c(-0.567, 0.048, 0.677),  # Correlation between intercept and slope for each class
    sd_r = 0.469,  # Residual standard deviation
    K = 3,  # Number of latent classes
    lambda_K = c(0.266, 0.269, 0.465)  # Class proportions
  )
  return(as.data.frame(UncGMM_dat))
}

# Load simulated datasets from specified directory
data_files <- list.files("~/Simulation_study/SimDat", full.names = TRUE, pattern = "SimulatedData_a")
SimDat <- data_files %>% map(read_csv)  # Load datasets into a list

# Function to load MCMC results from a specified directory
read_Stan_in_directory <- function(directory) {
  file_paths <- list.files(directory, full.names = TRUE, pattern = "Stan_a")
  file_contents <- map(file_paths, readRDS)
  return(file_contents)
}

# Function to load MCMC samples and diagnostic files from a specified subdirectory
read_files_in_directory <- function(subdirectory) {
  directory <- paste0("Simulation_study/MCMCResults_", subdirectory)
  
  # Load MCMC samples
  file_paths <- list.files(directory, full.names = TRUE, pattern = "rep_a")
  file_contents <- map(file_paths, readRDS)
  
  # Load log-likelihood and DI index files
  log_l <- readRDS(file.path(directory, "log_l.rds"))
  DI_data <- readRDS(file.path(directory, "DI_data.rds"))
  
  return(list(log_l = log_l, priors = file_contents, DI_data = DI_data))
}

# Load data and diagnostics for a specified example. Replace "D2C5" with your target subdirectory.
results_D2C5 <- read_files_in_directory("D2C5")
priors <- results_D2C5$priors
DI_data <- results_D2C5$DI_data

# Step 1: Rhat diagnostic on 4-chain batches
source("Diagnostics/Diagnostics.source.R")
Data_reordered_nonPermu <- traceData_ESS(priors, 1, ESS_var = "mu_intercept_1", ESS_chain = 50)$data  # Reorder data for ESS calculation
resulting_arrays <- split_data_into_arrays(Data_reordered_nonPermu, chains_per_array = 4)
Rhat_diag_by_chains(resulting_arrays, small_threshold = 1.05)  # Run Rhat diagnostic

# Generate traceplot for selected chains
num_chains <- 4  # Adjust to x if using the x-chain example dataset
iterations_per_chain <- 1000
total_iterations <- num_chains * iterations_per_chain
traceplot(Data_reordered = traceData(results_D2C5$priors, 1, total_iterations)$data, 
          num_chains = c(75:76), iterations_per_chain = 1000)$traceplot.by.chain  # Traceplot for chains 75 and 76

# Step 2: Stuck-sequence diagnostic to detect chains with persistent behavior
stuck <- lapply(seq_along(priors), 
                function(i) 
                  stuck_by_chain(priors[[i]], i, total_iter = total_iterations,
                                 iter_per_chain = 1000, window_size = 10, stuck_length = 20))
stuck  # Display stuck-sequence results

# Count chains exhibiting stuck behavior
stuck[[1]]$persistent_stuck_chain %>% length()
stuck[[1]]$num_chains_with_stuck

# Step 3: Twinlike-class diagnostic to detect similar class behavior across chains
twinlike_classes(DI_data = DI_data, 
                 selected_chains = 1:4, 
                 total_chains = num_chains, 
                 iter_per_chain = 1000, 
                 high_DI_value = 95, 
                 persist_length = 3,  
                 happen_times = 1)

# Step 4: Miniscule-class diagnostic to identify rare or small classes
diagnostics_graphs(
  Data = traceData(priors, 1, total_iterations),
  window_size = 10,
  selected_chains = 1:4,
  iter_per_chain = 1000,
  top_percentile_threshold = 0.95,
  DI_data = DI_data,
  total_chains = num_chains,
  high_DI_value = 95,
  persist_length = 3,
  happen_times = 1
)

