DI <- function(priors, num_classes, num_subjects, selected_chains, total_chains) {
  # Extract the mcmc draws
  #mcmc_draws <- priors[[1]][[2]][, selected_chains, -c(1:30)]
  mcmc_draws <- priors[[1]][[2]][, selected_chains, -c(1:21)]
  #mcmc_draws <- l[, selected_chains, ]
  # Calculate p_j^(k, l) for each subject, iteration, and chain
  num_iterations <- dim(mcmc_draws)[1]
  num_chains <- dim(mcmc_draws)[2]
  num_combinations <- num_classes * (num_classes - 1) / 2
  p_j_kl <- array(NA, dim = c(num_iterations, num_chains, num_subjects, num_combinations))
  entropy <- array(NA, dim = c(num_iterations, selected_chains%>%length(), num_subjects, num_combinations))
  # Initialize a vector to store column names
  column_names <- character(num_combinations)
  index <- 1
  
  for (iter in 1:num_iterations) {
    for (chain in 1:num_chains) {
      index <- 1
      for (j in 1:num_subjects) {
        for (k in 1:num_classes) {
          for (l in 1:num_classes) {
            if (k > l) {
              # Calculate p_j^(k, l) based on mcmc_draws
              logits <- c(mcmc_draws[iter, chain, (j-1)*num_classes + k], mcmc_draws[iter, chain, (j-1)*num_classes + l])
       
              # Apply softmax to get probabilities
              probabilities <- LDATS::softmax(logits)
              
              # Assign probabilities to the corresponding entries in p_j_kl
              p_j_kl[iter, chain, j, index] <- probabilities[1]  # Probability for class k
              p_j_kl[iter, chain, j, index] <- 
                ifelse(p_j_kl[iter, chain, j, index] > 0.5, 1 - p_j_kl[iter, chain, j, index], p_j_kl[iter, chain, j, index])
              
              
              entropy[iter, chain, j, index] <- probabilities[1] * log(probabilities[1] + 10e-99) +
                (1 - probabilities[1]) * log((1 - probabilities[1] + 10e-99))
              # Set column names
              column_names[index] <- sprintf("^(%d,%d)", k, l)
              index <- index + 1
            }
          }
        }
        index <- 1
      }
    }
  } 
  dimnames(p_j_kl) <- list(iteration = 1:num_iterations, chain = selected_chains,
                           subject = 1:num_subjects, combination = paste0("p_j", column_names))
  dimnames(entropy) <- list(iteration = 1:num_iterations, chain = selected_chains, 
                            subject = 1:num_subjects, combination = paste0("Entropy", column_names))
  
  
  average_p_kl <- apply(p_j_kl, c(1, 2, 4), mean)
  average_entropy <- - apply(entropy, c(1, 2, 4), mean)
  average_entropy3 <- -(average_p_kl * log(average_p_kl + 10e-99) +
                          (1 - average_p_kl) * log((1 - average_p_kl + 10e-99)))

  # Calculate distinguishability index for each pair (k, l) for each iteration and chain
  distinguishability_index <- (1 + (-1 / log(2)) * average_entropy) * 100
  distinguishability_index3 <- (1 + (-1 / log(2)) * average_entropy3) * 100
  x <- seq(1, total_chains * num_iterations)
  # Determine the starting and ending indices for each chain in num_chains
  indices <- lapply(selected_chains, function(chain) {
    start_index <- (chain - 1) * num_iterations + 1
    end_index <- chain * num_iterations
    return(start_index:end_index)
  })
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  x <- x[subset_indices]%>%as.data.frame()
  Data <- cbind(x, distinguishability_index[,,1]%>%array(, dim = c(num_iterations*num_chains, 1)),
                distinguishability_index[,,2]%>%array(, dim = c(num_iterations*num_chains, 1)),
                distinguishability_index[,,3]%>%array(, dim = c(num_iterations*num_chains, 1)))%>%as.data.frame()

  # Combine the indices and subset the data
  #subset_indices <- unlist(indices)
  #subset_data <- Data[subset_indices, ]%>%as.data.frame()
  
  # Generate permutations of length 2
  class_order <- order(as.data.frame(priors[[1]][[3]]$mean)[1:3,])
  permutations <- t(combn(class_order, 2))
  # Create column names based on permutations
  # colnames(Data) <- c("x",apply(permutations, 1, function(pair) paste0("DI", paste(pair, collapse = ""))))
  colnames(Data) <- c("x", 
                      paste(which(class_order == 1), "&", which(class_order == 2)), 
                      paste(which(class_order == 1), "&", which(class_order == 3)), 
                      paste(which(class_order == 2), "&", which(class_order == 3)))
  return(Data)
  
}
