This directory contains example datasets for running diagnostic functions, as demonstrated in ```Simulation_study/Sim_diagnostics_example.code.R```. These examples use 4 chains, whereas the analysis in the paper utilizes 100 chains. If using these example data files, please adjust the num_chains parameter in the diagnostics function accordingly, setting num_chains <- 4 for this example.

		rep_a1.rds: This file contains MCMC draws after post-hoc relabeling, generated using Simulation_study/Sim_design.code.R.

		log_l.rds: This file stores the log-likelihood values calculated post-relabeling (see Simulation_study/log_l.R for details on the calculation process). These files serve as 4-chain examples for the diagnostic function, but you can also generate them independently with your own configurations.

		DI_data.rds: This file contains the calculated Distinguishability Index (DI) values, introduced in the paper, using MCMC draws obtained after post-hoc relabeling (rep_a1.rds).

These example files allow for quick testing of diagnostic functions with smaller datasets and can be customized as needed.
