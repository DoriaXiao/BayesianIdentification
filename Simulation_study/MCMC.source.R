Stan_D6N50 = function(GMM_dat, K){
  GMMs_data_list_D6N50<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 1,
                     Normal_scale = 50,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6N50<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6N50(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T)
    #step_size = 0.01)
  return(flexSim_fit_3c_D6N50)
}


Stan_D6C5 = function(GMM_dat, K){
  GMMs_data_list_D6C5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 2,
                     Normal_scale = 5,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6C5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6C5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T)
    #step_size = 0.01)
  return(flexSim_fit_3c_D6C5)
}

Stan_D2C5_sss = function(GMM_dat, K){
  GMMs_data_list_D2C5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 2.0,
                     prior = 2,
                     Normal_scale = 5,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D2C5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D2C5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D2C5)
}
##############################################################################
Stan_D6N1_sss = function(GMM_dat, K){
  GMMs_data_list_D6N1<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 1,
                     Normal_scale = 1,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6N1<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6N1(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D6N1)
}


Stan_D6N5 = function(GMM_dat, K){
  GMMs_data_list_D6N5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 1,
                     Normal_scale = 5,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6N5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6N5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T)
    #step_size = 0.01)
  return(flexSim_fit_3c_D6N5)
}

Stan_D2N1_sss = function(GMM_dat, K){
  GMMs_data_list_D2N1<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 2.0,
                     prior = 1,
                     Normal_scale = 1,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D2N1<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D2N1(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D2N1)
}

Stan_D2N5 = function(GMM_dat, K){
  GMMs_data_list_D2N5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 2.0,
                     prior = 1,
                     Normal_scale = 5,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D2N5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D2N5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T)
    #step_size = 0.01)
  return(flexSim_fit_3c_D2N5)
}

#############################################################################
Stan_D10N50_sss = function(GMM_dat, K){
  GMMs_data_list_D10N50<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 10.0,
                     prior = 1,
                     Normal_scale = 50,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D10N50<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D10N50(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D10N50)
}



Stan_D6N50_sss = function(GMM_dat, K){
  GMMs_data_list_D6N50<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 1,
                     Normal_scale = 50,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6N50<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6N50(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D6N50)
}



Stan_D2N50_sss = function(GMM_dat, K){
  GMMs_data_list_D2N50<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 2.0,
                     prior = 1,
                     Normal_scale = 50,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D2N50<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D2N50(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D2N50)
}
##############################################################################
Stan_D10N100_sss = function(GMM_dat, K){
  GMMs_data_list_D10N100<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 10.0,
                     prior = 1,
                     Normal_scale = 100,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D10N100<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D10N100(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D10N100)
}



Stan_D6N100_sss = function(GMM_dat, K){
  GMMs_data_list_D6N100<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 1,
                     Normal_scale = 100,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6N100<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6N100(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D6N100)
}



Stan_D6C5_sss = function(GMM_dat, K){
  GMMs_data_list_D6C5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 6.0,
                     prior = 2,
                     Normal_scale = 100,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D6C5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D6C5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T,
    step_size = 0.01)
  return(flexSim_fit_3c_D6C5)
}

Stan_D2C5 = function(GMM_dat, K){
  GMMs_data_list_D2C5<- function(GMM_dat, K){
    cluster_size = GMM_dat %>% group_by(ID) %>% 
      dplyr::summarise(cluster_size = n()) %>% 
      pull(cluster_size)
    GMM_list <- list(Subject = as.numeric(factor(GMM_dat$ID, 
                                                 labels = 1:length(unique(GMM_dat$ID)))),
                     y=GMM_dat$y,
                     time=GMM_dat$time,
                     time_sq = GMM_dat$time_sq,
                     N = nrow(GMM_dat),
                     J = length(unique(GMM_dat$ID)),
                     s = cluster_size,
                     K=K,
                     Dir_alpha = 2.0,
                     prior = 2,
                     Normal_scale = 10,
                     Cauchy_scale = 5,
                     eta = 2
    )
    return(GMM_list)
  }
  flexSim_fit_3c_D2C5<- sq_GMM_ML_mod$sample(
    data = GMMs_data_list_D2C5(GMM_dat, K=K),
    chains = 100,
    parallel_chains = 4,
    iter_sampling = 1000,
    refresh = 1000,  
    save_warmup = T)
  #step_size = 0.01)
  return(flexSim_fit_3c_D2C5)
}



PP_sss = function(GMM_ML_fit_3c, chains, iterations){
  
  
  #GMM_ML_fit_3c<-flexSim_fit_3c_D6N3_sss
  K = 3
  m = chains*iterations
  S = 405
  J = 10
  
  var3<- c( "lambda" ,"mu_intercept", "mu_slope", "sq_mu_slope", 
            "sigma_e", "Omega[1,1,2]", "Omega[2,1,2]","Omega[3,1,2]", 
            "L_Omega[1,2,1]", "L_Omega[2,2,1]","L_Omega[3,2,1]",
            "L_Omega[1,2,2]", "L_Omega[2,2,2]","L_Omega[3,2,2]",
            "sigma_u[1,1]","sigma_u[1,2]","sigma_u[2,1]","sigma_u[2,2]","sigma_u[3,1]","sigma_u[3,2]", "pred_class_dis", "pred_class", "lp__")
  post_par_3c <- GMM_ML_fit_3c$draws(format = "df", variable = "pred_class") 
  post_par<- GMM_ML_fit_3c$draws(format = "df", variable = var3)
  
  post_par_3c_1 <- GMM_ML_fit_3c$draws(format = "df", variable = "pred_class")%>% 
    dplyr::select(ends_with('1]')) %>% as.matrix()
  post_par_3c_2 <- GMM_ML_fit_3c$draws(format = "df", variable = "pred_class")%>% 
    dplyr::select(ends_with('2]'))  %>% as.matrix()
  post_par_3c_3 <- GMM_ML_fit_3c$draws(format = "df", variable = "pred_class")%>% 
    dplyr::select(ends_with('3]'))  %>% as.matrix()
  post_class_p<-array (c (post_par_3c_1,post_par_3c_2, post_par_3c_3), dim=c (m,S,3))
  post_class <- array(data = NA, dim = c(m,S))
  for (i in 1:m){
    post_class[i,] <- apply(post_class_p[i,,], MARGIN = 1, which.max)
  }
  # initialize mcmc arrays
  mcmc <- array(data = NA, dim = c(m = m, K = K, J = J+S))
  
  
  mcmc[, 1, 1] <- post_par$`lambda[1]`
  mcmc[, 2, 1] <- post_par$`lambda[2]`
  mcmc[, 3, 1] <- post_par$`lambda[3]`
  mcmc[, 1, 2] <- post_par$`mu_intercept[1]`
  mcmc[, 2, 2] <- post_par$`mu_intercept[2]`
  mcmc[, 3, 2] <- post_par$`mu_intercept[3]`
  mcmc[, 1, 3] <- post_par$`mu_slope[1]`
  mcmc[, 2, 3] <- post_par$`mu_slope[2]`
  mcmc[, 3, 3] <- post_par$`mu_slope[3]`
  mcmc[, 1, 4] <- post_par$`sq_mu_slope[1]`
  mcmc[, 2, 4] <- post_par$`sq_mu_slope[2]`
  mcmc[, 3, 4] <- post_par$`sq_mu_slope[3]`
  mcmc[, 1, 5] <- post_par$`Omega[1,1,2]`
  mcmc[, 2, 5] <- post_par$`Omega[2,1,2]`
  mcmc[, 3, 5] <- post_par$`Omega[3,1,2]`
  mcmc[, 1, 6] <- post_par$`sigma_u[1,1]`
  mcmc[, 2, 6] <- post_par$`sigma_u[2,1]`
  mcmc[, 3, 6] <- post_par$`sigma_u[3,1]`
  mcmc[, 1, 7] <- post_par$`sigma_u[1,2]`
  mcmc[, 2, 7] <- post_par$`sigma_u[2,2]`
  mcmc[, 3, 7] <- post_par$`sigma_u[3,2]`
  mcmc[, , 8] <- post_par$sigma_e
  mcmc[, 1, 9] <- post_par$`L_Omega[1,2,1]`
  mcmc[, 2, 9] <- post_par$`L_Omega[2,2,1]`
  mcmc[, 3, 9] <- post_par$`L_Omega[3,2,1]`
  mcmc[, 1, 10] <- post_par$`L_Omega[1,2,2]`
  mcmc[, 2, 10] <- post_par$`L_Omega[2,2,2]`
  mcmc[, 3, 10] <- post_par$`L_Omega[3,2,2]`
  for (i in 1:S){
    mcmc[,1,i+10] <-as.matrix(post_par)[, paste('pred_class[',i,',1]', sep = '')]
    mcmc[,2,i+10] <-as.matrix(post_par)[, paste('pred_class[',i,',2]', sep = '')]
    mcmc[,3,i+10] <-as.matrix(post_par)[, paste('pred_class[',i,',3]', sep = '')]
  }
  
  
  source("PostProcessing_list.R")
  fit_permuted_3c <-post_processing(chains = chains, iterations = iterations, K=3, J=J+S, post_class, mcmc, post_class_p, post_par)
  return(fit_permuted_3c)
}
