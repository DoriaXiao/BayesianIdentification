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

