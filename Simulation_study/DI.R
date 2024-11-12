DI <- function(priors, num_classes, num_subjects, selected_chains, total_chains) {
  # Extract the MCMC draws of parameters of interest
  mcmc_draws <- priors[[1]][[2]][, selected_chains, -c(1:21)]
  
  # Initialize arrays for probabilities (p_j^(k, l)) and entropy
  num_iterations <- dim(mcmc_draws)[1]
  num_chains <- dim(mcmc_draws)[2]
  num_combinations <- num_classes * (num_classes - 1) / 2
  p_j_kl <- array(NA, dim = c(num_iterations, num_chains, num_subjects, num_combinations))
  entropy <- array(NA, dim = c(num_iterations, length(selected_chains), num_subjects, num_combinations))
  
  # Initialize column names and index
  column_names <- character(num_combinations)
  
  for (iter in 1:num_iterations) {
    for (chain in 1:num_chains) {
      index <- 1
      for (j in 1:num_subjects) {
        for (k in 1:num_classes) {
          for (l in 1:num_classes) {
            if (k > l) {
              # Calculate logits and apply softmax to obtain probabilities
              logits <- c(mcmc_draws[iter, chain, (j - 1) * num_classes + k], 
                          mcmc_draws[iter, chain, (j - 1) * num_classes + l])
              probabilities <- LDATS::softmax(logits)
              
              # Calculate and adjust p_j^(k, l) to lie between 0 and 0.5
              p_j_kl[iter, chain, j, index] <- ifelse(probabilities[1] > 0.5, 
                                                      1 - probabilities[1], 
                                                      probabilities[1])
              
              # Calculate entropy for p_j^(k, l)
              entropy[iter, chain, j, index] <- p_j_kl[iter, chain, j, index] * log(p_j_kl[iter, chain, j, index] + 1e-99) +
                                                (1 - p_j_kl[iter, chain, j, index]) * log(1 - p_j_kl[iter, chain, j, index] + 1e-99)
              
              # Set column names for p_j^(k, l) combinations
              column_names[index] <- sprintf("^(%d,%d)", k, l)
              index <- index + 1
            }
          }
        }
      }
    }
  }
  
  # Add dimension names to arrays
  dimnames(p_j_kl) <- list(iteration = 1:num_iterations, chain = selected_chains,
                           subject = 1:num_subjects, combination = paste0("p_j", column_names))
  dimnames(entropy) <- list(iteration = 1:num_iterations, chain = selected_chains, 
                            subject = 1:num_subjects, combination = paste0("Entropy", column_names))
  
  # Calculate average probabilities and entropy
  average_p_kl <- apply(p_j_kl, c(1, 2, 4), mean)
  average_entropy <- -apply(entropy, c(1, 2, 4), mean)
  average_entropy3 <- -(average_p_kl * log(average_p_kl + 1e-99) +
                          (1 - average_p_kl) * log(1 - average_p_kl + 1e-99))
  
  # Calculate the distinguishability index
  distinguishability_index <- (1 + (-1 / log(2)) * average_entropy) * 100
  distinguishability_index3 <- (1 + (-1 / log(2)) * average_entropy3) * 100
  
  # Create an index sequence for chains
  x <- seq(1, total_chains * num_iterations)
  indices <- lapply(selected_chains, function(chain) {
    start_index <- (chain - 1) * num_iterations + 1
    end_index <- chain * num_iterations
    start_index:end_index
  })
  subset_indices <- unlist(indices)
  x <- as.data.frame(x[subset_indices])
  
  # Create data frame for the distinguishability index
  Data <- cbind(x, 
                distinguishability_index[,,1] %>% array(dim = c(num_iterations * num_chains, 1)),
                distinguishability_index[,,2] %>% array(dim = c(num_iterations * num_chains, 1)),
                distinguishability_index[,,3] %>% array(dim = c(num_iterations * num_chains, 1))) %>% 
    as.data.frame()
  
  # Generate class order and create column names based on class combinations
  class_order <- order(as.data.frame(priors[[1]][[3]]$mean)[1:3,])
  colnames(Data) <- c("x", 
                      paste(which(class_order == 1), "&", which(class_order == 2)), 
                      paste(which(class_order == 1), "&", which(class_order == 3)), 
                      paste(which(class_order == 2), "&", which(class_order == 3)))
  
  return(Data)
}
