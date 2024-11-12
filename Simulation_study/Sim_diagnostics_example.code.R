# Simulate data
setwd("C:/Users/xingy/OneDrive/2023 Fall/BayesIdentify/GitHub_code")
source("~/Simulation_study/SimCode.source.R")

data_fun_MCMC <- function() {
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
  return(as.data.frame(UncGMM_dat))
}

# Load saved simulated datasets
data_files <- list.files("~/Simulation_study/SimDat", full.names = TRUE, pattern = "SimulatedData_a")
SimDat <- data_files %>% map(read_csv)

# Functions to load MCMC results
read_Stan_in_directory <- function(directory) {
  file_paths <- list.files(directory, full.names = TRUE, pattern = "Stan_a")
  file_contents <- map(file_paths, readRDS)
  return(file_contents)
}

read_files_in_directory <- function(subdirectory) {
  directory <- paste("Simulation_study/MCMCResults_", subdirectory, sep = "")
  file_paths <- list.files(directory, full.names = TRUE, pattern = "rep_a")
  file_contents <- map(file_paths, readRDS)
  
  # Read the log_l.rds file
  log_l_path <- file.path(directory, "log_l.rds")
  log_l <- readRDS(log_l_path)
  
  return(list(log_l = log_l, priors = file_contents))
}

# Replace "D2C5" with your desired subdirectory
results_D2C5 <- read_files_in_directory("D2C5")

# Load DI index results
DI_results_D2 <- readRDS("Simulation_study/Sim_result_example/DI_results_D2.rds")

priors <- results_D2C5$priors

# Step 1: Rhat diagnostic by 4-chain batches
source("Diagnostics/Diagnostics.source.R")
Data_reordered_nonPermu <- traceData_ESS(priors, 1, ESS_var = "mu_intercept_1", ESS_chain = 50)$data # [iterations, chains, parameters]
resulting_arrays <- split_data_into_arrays(Data_reordered_nonPermu, chains_per_array = 4)
Rhat_diag_by_chains(resulting_arrays, small_threshold = 1.05)


# Traceplot by chains
num_chains <- 100
iterations_per_chain <- 1000
total_iterations <- num_chains*iterations_per_chain
traceplot(Data_reordered = traceData(results_D2C5$priors, 1, total_iterations)$data, 
          num_chains = c(75:76), iterations_per_chain = 1000)$traceplot.by.chain

# Step 2: Stuck-sequence diagnostic
stuck <- lapply(seq_along(priors), 
                function(i) 
                  stuck_by_chain(priors[[i]], i, total_iter = total_iterations,
                                 iter_per_chain = 1000, window_size = 10, stuck_length = 20))
stuck

stuck[[1]]$persistent_stuck_chain%>%length()
stuck[[1]]$num_chains_with_stuck

# Step 3: Twinlike-class diagnostic
twinlike_classes(DI_data = DI_data, 
                 selected_chains = c(1:100), 
                 total_chains = num_chains, 
                 iter_per_chain = 1000, 
                 high_DI_value = 95, 
                 persist_length = 3,  
                 happen_times = 1)

# Step 4: Miniscule-class diagnostic example
DI_data = DI_results_D2$D2C5
diagnostics_graphs(
  Data = traceData(priors, 1, total_iterations),
  window_size = 10,
  selected_chains = c(1:100),
  iter_per_chain = 1000,
  top_percentile_threshold = 0.95,
  DI_data = DI_data,
  total_chains = num_chains,
  high_DI_value = 95,
  persist_length = 3,
  happen_times = 1
)


