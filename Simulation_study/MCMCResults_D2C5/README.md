# Example Datasets for Diagnostics

This directory contains example datasets for running diagnostic functions, as demonstrated in [Sim_diagnostics_example.code.R](../../Sim_diagnostics_example.code.R). These examples use 4 chains, whereas the analysis in the paper utilizes 100 chains. If using these example data files, please adjust the `num_chains` parameter in the diagnostics function accordingly, setting `num_chains <- 4` for this example.

- **`rep_a1.rds`**: Contains MCMC draws after post-hoc relabeling, generated with [Sim_design.code.R](../../Sim_design.code.R).

- **`log_l.rds`**: Stores log-likelihood values calculated after relabeling (refer to [log_l.R](../../Simulation_study/log_l.R) for the computation process). This file serves as a 4-chain example for the diagnostic function, but can also be generated independently with different configurations.

- **`DI_data.rds`**: Contains the calculated Distinguishability Index (DI) values, as introduced in the paper, using MCMC draws from `rep_a1.rds` after post-hoc relabeling.

These example files allow for quick testing of diagnostic functions with smaller datasets and are customizable as needed.
