# R set up
library(haven)
library(gridExtra)
library(magrittr)
library(tidyverse) # for ggplot2, dplyr, purrr, etc.
library(MASS) # for mvrnorm
library(fs) 
library("loo")
# Installing CmdStan
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
library(tictoc)
library(furrr)
plan("multisession")
library(PropCIs) # for calculating confidence intervals for a binomial proportion
color_scheme_set("brightblue")
library(kableExtra)
#install_cmdstan(cores = 16)
cmdstan_path()
cmdstan_version()
library(flexmix)
library(rstan)
setwd("~/Simulation_study")

#################
##D2##6##10######
##C5##N100##N50##
##D##SS##########
#################


source("~/SimCode.source.R")
data_fun_MCMC<-function(){
  UncGMM_dat <- UncGMM_data(n_pers=405,
                            n_time=4,
                            beta_int=c(3.417, 2.637, 2.019), 
                            beta_slo_time=c(1.297, 2.607, 1.393), 
                            beta_slo_time_sq=c(-0.094, -0.472, -0.120),
                            sd_i=c(0.967, 0.599, 0.188), 
                            sd_s=c(0.263, 0.265, 0.356),  
                            cor_is=c(-0.567, 0.048, 0.677), 
                            sd_r=0.469, 
                            K = 3,
                            lambda_K = c(0.266, 0.269, 0.465))
  return(as.data.frame(UncGMM_dat))
}


# save simulated datasets (this has been saved to "~/SimDat")
#sapply(1:1, function(i) write.csv(data_MCMC[i], sprintf("SimulatedData_a%i.csv",i)))

# load saved simulated datasets

data_files<- list.files("~/SimDat", 
                        full.names = TRUE,
                        pattern="SimulatedData_a",)  # Identify file names
SimDat<- data_files%>% 
  map(function(path){
    read_csv(path)
  })   

# plot simulated data
ggplot(SimDat[[1]], aes(x = time, y = y, 
                        group = ID%>%as.factor(), 
                        color = latent_group%>%as.factor())) +
  geom_line() +
  scale_color_manual(values = c("red", "blue", "black"))  


# 2. Running model
# Compiling a GMM with marginal likelihood
sq_GMM_ML_file <- file.path(cmdstan_path(), "GMM_ML.stan")
sq_GMM_ML_mod <- cmdstan_model(sq_GMM_ML_file)

###############################D2C5#################################################

for (i in 1:1) {
  source("~/MCMC.source.R")
  setwd("~/BayesIdentify/MCMCResults_D2C5")
  Stan_D2C5(SimDat[[i]], K=3)$save_object( sprintf("Stan_a%i.rds",i))
}

# load MCMC results

Stan_paths_rds_D2C5 <- list.files("~/MCMCResults_D2C5", 
                                     full.names = TRUE,
                                     pattern="Stan_a",)  # Identify file names
Stan_contents_rds_D2C5 <- Stan_paths_rds_D2C5%>% 
  map(function(path){
    readRDS(path)
  })


# run replications

for (i in 1:1) {
  source("~/MCMC.source.R")
  setwd("~/MCMCResults_D2C5")
  saveRDS(PP_sss(Stan_contents_rds_D2C5[[i]], 100, 1000), sprintf("rep_a%i.rds",i))
}

# load pp results

file_paths_rds_D2C5 <- list.files("~/MCMCResults_D2C5", 
                                     full.names = TRUE,
                                     pattern="rep_a",)  # Identify file names
file_contents_rds_D2C5 <- file_paths_rds_D2C5%>% 
  map(function(path){
    readRDS(path)
  })


