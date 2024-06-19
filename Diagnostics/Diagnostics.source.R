library(RcppRoll)
library(ggplot2)
library(tidyr)
library(patchwork)
library(cowplot)
library(LDATS)
library(coda)
library(abind)
library(mixtools)
library(zoo)
library(grid)
traceData_ESS <- function(file_contents_rds_, i, ESS_var, ESS_chain){
  lambda_1 <- file_contents_rds_[[i]][[2]][,,1]
  sigma1_1 <- file_contents_rds_[[i]][[2]][,,13]
  sigma2_1 <- file_contents_rds_[[i]][[2]][,,16]
  lambda_2 <- file_contents_rds_[[i]][[2]][,,2]
  sigma1_2 <- file_contents_rds_[[i]][[2]][,,14]
  sigma2_2 <- file_contents_rds_[[i]][[2]][,,17]
  lambda_3 <- file_contents_rds_[[i]][[2]][,,3]
  sigma1_3 <- file_contents_rds_[[i]][[2]][,,15]
  sigma2_3 <- file_contents_rds_[[i]][[2]][,,18]
  Omega12_1 <- file_contents_rds_[[i]][[2]][,,10]
  Omega12_2 <- file_contents_rds_[[i]][[2]][,,11]
  Omega12_3 <- file_contents_rds_[[i]][[2]][,,12]
  mu_intercept_1<-file_contents_rds_[[i]][[2]][,,4]
  mu_intercept_2<-file_contents_rds_[[i]][[2]][,,5]
  mu_intercept_3<-file_contents_rds_[[i]][[2]][,,6]
  mu_slope_1<-file_contents_rds_[[i]][[2]][,,7]
  mu_slope_2<-file_contents_rds_[[i]][[2]][,,8]
  mu_slope_3<-file_contents_rds_[[i]][[2]][,,9]
  sq_mu_slope_1<-file_contents_rds_[[i]][[2]][,,22]
  sq_mu_slope_2<-file_contents_rds_[[i]][[2]][,,23]
  sq_mu_slope_3<-file_contents_rds_[[i]][[2]][,,24]
  # Assuming lambda_1, lambda_2, and lambda_3 are the means of lambda for each class
  class_order <- order(c(mean(lambda_1), mean(lambda_2), mean(lambda_3)))
  
  # Define the order of variables for each class
  lambda_cols <- c("lambda", "sigma1", "sigma2")
  mu_intercept_cols <- c("mu_intercept", "mu_slope", "sq_mu_slope")
  omega_cols <- c("Omega12")
  
  
  # Create the variable names for each class based on the order
  lambda_var_names <- c(sapply(lambda_cols, function(col) paste(col, class_order, sep = "_")))
  mu_intercept_var_names <- c(sapply(mu_intercept_cols, function(col) paste(col, class_order, sep = "_")))
  omega_var_names <- c(sapply(omega_cols, function(col) paste(col, class_order, sep = "_")))
  
  # Combine all variable names
  all_var_names <- c(lambda_var_names, mu_intercept_var_names, omega_var_names)
  
  Data <- abind(lambda_1 = lambda_1, sigma1_1 = sigma1_1, sigma2_1= sigma2_1,
                lambda_2 = lambda_2, sigma1_2 = sigma1_2, sigma2_2= sigma2_2,
                lambda_3 = lambda_3, sigma1_3 = sigma1_3, sigma2_3= sigma2_3,
                mu_intercept_1 = mu_intercept_1, mu_intercept_2 = mu_intercept_2, mu_intercept_3 = mu_intercept_3,
                mu_slope_1 = mu_slope_1, mu_slope_2 = mu_slope_2, mu_slope_3 = mu_slope_3,
                sq_mu_slope_1 = sq_mu_slope_1, sq_mu_slope_2 = sq_mu_slope_2, sq_mu_slope_3 = sq_mu_slope_3,
                Omega12_1 = Omega12_1, Omega12_2 = Omega12_2, Omega12_3 = Omega12_3,  along = 3)
  
  return_list<-list(
    Data_reordered <- abind(lambda_1 = Data[,,all_var_names[1]], sigma1_1 = Data[,,all_var_names[4]], sigma2_1= Data[,,all_var_names[7]],
                            lambda_2 = Data[,,all_var_names[2]], sigma1_2 = Data[,,all_var_names[5]], sigma2_2= Data[,,all_var_names[8]],
                            lambda_3 = Data[,,all_var_names[3]], sigma1_3 = Data[,,all_var_names[6]], sigma2_3= Data[,,all_var_names[9]],
                            mu_intercept_1 = Data[,,all_var_names[10]], mu_intercept_2 = Data[,,all_var_names[11]], mu_intercept_3 = Data[,,all_var_names[12]],
                            mu_slope_1 = Data[,,all_var_names[13]], mu_slope_2 = Data[,,all_var_names[14]], mu_slope_3 = Data[,,all_var_names[15]],
                            sq_mu_slope_1 = Data[,,all_var_names[16]], sq_mu_slope_2 = Data[,,all_var_names[17]], sq_mu_slope_3 = Data[,,all_var_names[18]],
                            Omega12_1 = Data[,,all_var_names[19]], Omega12_2 = Data[,,all_var_names[20]], Omega12_3 = Data[,,all_var_names[21]], along = 3)
  )
  return_list$data <- Data_reordered
  file_contents <-Data_reordered[,,ESS_var]
  rownames(file_contents) <- NULL
  colnames(file_contents) <- NULL
  create_mcmc_list <- function(mcmc_matrix) {
    num_chains <- dim(mcmc_matrix)[2]
    mcmc_chain_list <- vector("list", length = num_chains)
    
    for (i in 1:num_chains) {
      chain_i <- mcmc_matrix[, i, drop = FALSE]
      mcmc_chain_list[[i]] <- mcmc(chain_i)
    }
    
    return(mcmc.list(mcmc_chain_list))
  }
  mcmc_list <- create_mcmc_list(file_contents)
  return_list$coda.ESS.chain<-coda::effectiveSize(mcmc_list[[ESS_chain]])
  return_list$coda.ESS.total<-coda::effectiveSize(mcmc_list)
  return_list$coda.autocorr.chain<-coda::autocorr(mcmc_list[[ESS_chain]])
  return_list$coda.ESS.allchain <-lapply(seq_along(mcmc_list), 
                                         function(i) effectiveSize(mcmc_list[[i]]))%>%as.matrix()
  return(return_list)
}


traceData<-function(file_contents_rds_, i, iter){
  lambda_1 <- file_contents_rds_[[i]][[2]][,,1]%>%array(, dim = c(iter, 1))
  sigma1_1 <- file_contents_rds_[[i]][[2]][,,13]%>%array(, dim = c(iter, 1))
  sigma2_1 <- file_contents_rds_[[i]][[2]][,,16]%>%array(, dim = c(iter, 1))
  lambda_2 <- file_contents_rds_[[i]][[2]][,,2]%>%array(, dim = c(iter, 1))
  sigma1_2 <- file_contents_rds_[[i]][[2]][,,14]%>%array(, dim = c(iter, 1))
  sigma2_2 <- file_contents_rds_[[i]][[2]][,,17]%>%array(, dim = c(iter, 1))
  lambda_3 <- file_contents_rds_[[i]][[2]][,,3]%>%array(, dim = c(iter, 1))
  sigma1_3 <- file_contents_rds_[[i]][[2]][,,15]%>%array(, dim = c(iter, 1))
  sigma2_3 <- file_contents_rds_[[i]][[2]][,,18]%>%array(, dim = c(iter, 1))
  Omega12_1 <- file_contents_rds_[[i]][[2]][,,10]%>%array(, dim = c(iter, 1))
  Omega12_2 <- file_contents_rds_[[i]][[2]][,,11]%>%array(, dim = c(iter, 1))
  Omega12_3 <- file_contents_rds_[[i]][[2]][,,12]%>%array(, dim = c(iter, 1))
  mu_intercept_1<-file_contents_rds_[[i]][[2]][,,4]%>%array(, dim = c(iter, 1))
  mu_intercept_2<-file_contents_rds_[[i]][[2]][,,5]%>%array(, dim = c(iter, 1))
  mu_intercept_3<-file_contents_rds_[[i]][[2]][,,6]%>%array(, dim = c(iter, 1))
  mu_slope_1<-file_contents_rds_[[i]][[2]][,,7]%>%array(, dim = c(iter, 1))
  mu_slope_2<-file_contents_rds_[[i]][[2]][,,8]%>%array(, dim = c(iter, 1))
  mu_slope_3<-file_contents_rds_[[i]][[2]][,,9]%>%array(, dim = c(iter, 1))
  sq_mu_slope_1<-file_contents_rds_[[i]][[2]][,,22]%>%array(, dim = c(iter, 1))
  sq_mu_slope_2<-file_contents_rds_[[i]][[2]][,,23]%>%array(, dim = c(iter, 1))
  sq_mu_slope_3<-file_contents_rds_[[i]][[2]][,,24]%>%array(, dim = c(iter, 1))
  # Assuming lambda_1, lambda_2, and lambda_3 are the means of lambda for each class
  class_order <- order(c(mean(lambda_1), mean(lambda_2), mean(lambda_3)))
  
  # Define the order of variables for each class
  lambda_cols <- c("lambda", "sigma1", "sigma2")
  mu_intercept_cols <- c("mu_intercept", "mu_slope", "sq_mu_slope")
  omega_cols <- c("Omega12")
  
  
  # Create the variable names for each class based on the order
  lambda_var_names <- c(sapply(lambda_cols, function(col) paste(col, class_order, sep = "_")))
  mu_intercept_var_names <- c(sapply(mu_intercept_cols, function(col) paste(col, class_order, sep = "_")))
  omega_var_names <- c(sapply(omega_cols, function(col) paste(col, class_order, sep = "_")))
  
  # Combine all variable names
  all_var_names <- c(lambda_var_names, mu_intercept_var_names, omega_var_names)
  
  x <- seq(1, iter)
  Data <- data.frame(x = x, lambda_1 = lambda_1, sigma1_1 = sigma1_1, sigma2_1= sigma2_1,
                     lambda_2 = lambda_2, sigma1_2 = sigma1_2, sigma2_2= sigma2_2,
                     lambda_3 = lambda_3, sigma1_3 = sigma1_3, sigma2_3= sigma2_3,
                     mu_intercept_1 = mu_intercept_1, mu_intercept_2 = mu_intercept_2, mu_intercept_3 = mu_intercept_3,
                     mu_slope_1 = mu_slope_1, mu_slope_2 = mu_slope_2, mu_slope_3 = mu_slope_3,
                     sq_mu_slope_1 = sq_mu_slope_1, sq_mu_slope_2 = sq_mu_slope_2, sq_mu_slope_3 = sq_mu_slope_3,
                     Omega12_1 = Omega12_1, Omega12_2 = Omega12_2, Omega12_3 = Omega12_3)
  
  return_list<-list(
    Data_reordered <- data.frame(x = x, lambda_1 = Data[[all_var_names[1]]], sigma1_1 = Data[[all_var_names[4]]], sigma2_1= Data[[all_var_names[7]]],
                                 lambda_2 = Data[[all_var_names[2]]], sigma1_2 = Data[[all_var_names[5]]], sigma2_2= Data[[all_var_names[8]]],
                                 lambda_3 = Data[[all_var_names[3]]], sigma1_3 = Data[[all_var_names[6]]], sigma2_3= Data[[all_var_names[9]]],
                                 mu_intercept_1 = Data[[all_var_names[10]]], mu_intercept_2 = Data[[all_var_names[11]]], mu_intercept_3 = Data[[all_var_names[12]]],
                                 mu_slope_1 = Data[[all_var_names[13]]], mu_slope_2 = Data[[all_var_names[14]]], mu_slope_3 = Data[[all_var_names[15]]],
                                 sq_mu_slope_1 = Data[[all_var_names[16]]], sq_mu_slope_2 = Data[[all_var_names[17]]], sq_mu_slope_3 = Data[[all_var_names[18]]],
                                 Omega12_1 = Data[[all_var_names[19]]], Omega12_2 = Data[[all_var_names[20]]], Omega12_3 = Data[[all_var_names[21]]])
  )
  return_list$data <- Data_reordered
  return_list$traceplot<- ggplot(Data_reordered, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", size = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", size = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", size = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( expression(paste(lambda^(3))), 
                                   expression(paste(lambda^(2))), 
                                   expression(paste(lambda^(1)))))
  
  return(return_list)
}




find_positions <- function(x, small_threshold, large_threshold) {
  # Find the positions of elements that are smaller than the small threshold
  x<- x$diagnostic_summary()$ebfmi
  small_positions <- which(x < small_threshold)
  # Find the positions of elements that are larger than the large threshold
  large_positions <- which(x > large_threshold)
  
  # Print a warning message if there are any small positions
  if (length(small_positions) > 0) {
    warning(paste("The following chains had an E-BFMI less than 0.3:", 
                  paste(small_positions, collapse = ", ")))
  } else {
    # Return a message if there are no small positions
    print("No chains had an E-BFMI less than 0.3.")
  }
  
  # Print a warning message if there are any large positions
  if (length(large_positions) > 0) {
    warning(paste("The following chains had an E-BFMI greater than 1.3:", 
                  paste(large_positions, collapse = ", ")))
  } else {
    # Return a message if there are no large positions
    print("No chains had an E-BFMI greater than 1.3.")
  }
  
  # Return the small positions
  #return(small_positions)
}

find_positions_cross <- function(x, small_threshold, large_threshold, index) {
  # Find the positions of elements that are smaller than the small threshold
  x<-x$diagnostic_summary()$ebfmi
  small_positions <- which(x < small_threshold)
  # Find the positions of elements that are larger than the large threshold
  large_positions <- which(x > large_threshold)
  
  # Initialize an empty list to store the warning messages and positions
  warnings <- list()
  positions <- list()
  
  # Store the warning messages and positions in the lists
  if (length(small_positions) > 0) {
    warnings$small <- paste("The following chains had an E-BFMI less than", small_threshold, ":",
                            paste(small_positions, collapse = ", "))
    positions$small <- small_positions
  }
  if (length(large_positions) > 0) {
    warnings$large <- paste("The following chains had an E-BFMI greater than", large_threshold, ":",
                            paste(large_positions, collapse = ", "))
    positions$large <- large_positions
  }
  
  # Return the list of warning messages and positions
  return(list(warnings = Filter(Negate(is.null), warnings), index = index))
}

# This function identifies the positions of Rhat values greater than a small threshold
find_rhat_mean_cross <- function(x, small_threshold, index) {
  
  # Compute the mean Rhat value across all parameters
  rhat_mean <- mean(summarise_draws(x[[2]][,,1:21])$rhat)
  
  # Check if the Rhat mean exceeds the small threshold
  if (rhat_mean > small_threshold) {
    
    # Generate a warning message for Rhat values that exceed the threshold
    warning_msg <- paste("The Rhat values for the parameters in your analysis are higher than the recommended threshold of", small_threshold, 
                         "with a mean of", round(rhat_mean, 2))
    
    # Return the warning message and the input index
    return(list(warnings = list(large_rhat = warning_msg), index = index))
    
  } else {
    
    # Generate a message indicating that the Rhat values are within the threshold
    message <- paste("The Rhat values for the parameters in your analysis are not higher than the recommended threshold of", small_threshold)
    
    # Return the message and the input index
    return(list(message = message, index = index))
  }
}




# This function identifies the positions of Rhat values greater than a small threshold
find_rhat_positions_cross <- function(x, small_threshold, index) {
  position <- c(summarise_draws(x[[2]][,,1:19])$rhat , summarise_draws(x[[2]][,,22:24])$rhat)
  # Identify the positions with Rhat values greater than the threshold
  large_rhat_positions <- which(position > small_threshold)
  larger_rhat_positions <- which(position > 1.1)
  
  # Initialize an empty list to store the warning messages and positions
  warnings <- list()
  positions <- list()
  
  # If there are positions with large Rhat values
  if (length(large_rhat_positions) > 0) {
    
    # Create a warning message
    warnings$large_rhat <- paste("The following parameter positions had a Rhat greater than", small_threshold, ":", 
                                 paste(large_rhat_positions, collapse = ", "))
    
    # Store the positions in the list
    positions$large_rhat <- large_rhat_positions
    
  } else {
    # If there are no positions with large Rhat values, print a message
    print(paste("Simulation replicate ", index ," has NO parameters that had a Rhat greater than ", small_threshold, "."))
  }
  
  # Print the positions in the list
  if (length(positions$large_rhat) > 0) {
    
    message("Simulation replicate ", index, " exhibits ", (large_rhat_positions %>% length() - 1), " out of 21 parameters with Rhat > ", small_threshold, 
            ". The average Rhat across all 21 parameters is ", position %>% mean, ", and there are ", (larger_rhat_positions %>% length() - 1), " parameters with Rhat > 1.10.")
    
    
  }
  
  # Return the list of warnings and the index
  return(list(warnings = Filter(Negate(is.null), warnings), positions = Filter(Negate(is.null), positions), index = index))
}



# 3. % of time lambdas less than 1%

small_lambda <-function(file_contents_rds_, i){
  sum<-sum(file_contents_rds_[[i]][[2]][,,1] < 0.01 |
             file_contents_rds_[[i]][[2]][,,2] < 0.01 |
             file_contents_rds_[[i]][[2]][,,3] < 0.01) /
    (file_contents_rds_[[i]][[2]][,,1]%>%length()*3)
  sum_percent <- round(sum*100, 2)
  
  if (sum > 0.0001) {
    warning(paste("Outlier class draws occur in",paste(sum_percent, collapse = ""),
                  "% of the iterations."))
  } else {
    # Return a message if there are no small positions
    return("No outlier class was detected in any of the iterations.")
  }
  #return(sum)
}

determine_msd_threshold <- function(Data, window_size){
  mean_lambda_1 <- mean(Data$lambda_1)
  mean_lambda_2 <- mean(Data$lambda_2)
  mean_lambda_3 <- mean(Data$lambda_3)
  
  # Determine which variable has the smallest mean
  if (mean_lambda_1 < mean_lambda_2 && mean_lambda_1 < mean_lambda_3) {
    smallest_mean <- "lambda_1"
  } else if (mean_lambda_2 < mean_lambda_1 && mean_lambda_2 < mean_lambda_3) {
    smallest_mean <- "lambda_2"
  } else {
    smallest_mean <- "lambda_3"
  }
  
  Data$smallest_mean<-Data[[smallest_mean]]

  mov_sd <- RcppRoll::roll_sd(Data$smallest_mean, window_size, fill = NA)

  fit<-normalmixEM(na.omit(mov_sd), k = 2)
  return(list(summary <- cbind(lambda <- fit$lambda, mu <- fit$mu, sigma <- fit$sigma), 
              threshold <- fit$mu[1]))
}

diagnostic_moving_sd<- function(Data,  window_size, change_threshold,tiny_value_threshold) {
  
  # Calculate the mean of each variable
  mean_lambda_1 <- mean(Data$lambda_1)
  mean_lambda_2 <- mean(Data$lambda_2)
  mean_lambda_3 <- mean(Data$lambda_3)
  
  # Determine which variable has the smallest mean
  if (mean_lambda_1 < mean_lambda_2 && mean_lambda_1 < mean_lambda_3) {
    smallest_mean <- "lambda_1"
  } else if (mean_lambda_2 < mean_lambda_1 && mean_lambda_2 < mean_lambda_3) {
    smallest_mean <- "lambda_2"
  } else {
    smallest_mean <- "lambda_3"
  }
  
  Data$smallest_mean<-Data[[smallest_mean]]

  mov_sd <- RcppRoll::roll_sd(Data$smallest_mean, window_size, fill = NA)
  
  df <- data.frame(x = seq_along(mov_sd), y = mov_sd)
  # Plot the moving standard deviation using ggplot
  p<-ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    labs(x = "Iteration", y = "Moving Standard Deviation", 
         title = "")+
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray", linetype = "dotted"),
      panel.grid.minor = element_blank())
  
  
  # Calculate the maximum difference between consecutive elements in the vector
  max_diff <- max(abs(diff(mov_sd)), na.rm = TRUE)
  min <- min(mov_sd, na.rm = TRUE)
  message <- paste("The variable with the smallest mean is", smallest_mean, ".\n")
  if (max_diff > change_threshold) {
    if (min > tiny_value_threshold) {
      message <- paste(message, "There is a sudden change in moving standard deviation:", max_diff, ".\n")
    } else {
      message <- paste(message, "There is a sudden change in moving standard deviation:", max_diff, ".\n")
      message <- paste(message, "The moving standard deviation stays at a tiny value:", min, ".\n")
      #message <- paste(message, "The moving standard deviation stays at zero for ",  which(abs(diff(mov_sd))==0)%>%length(), " iterations", ".\n")
    }
  } else {
    if (min > tiny_value_threshold) {
      message <- paste(message, "No sudden change in moving standard deviation. The maximum change is:", max_diff, ".\n")
    } else {
      message <- paste(message, "No sudden change in moving standard deviation. The maximum change is:", max_diff, ".\n")
      message <- paste(message, "The moving standard deviation stays at a tiny value:", min, ".\n")
      #message <- paste(message, "The moving standard deviation stays at zero for ",  which(abs(diff(mov_sd))==0)%>%length(), " iterations", ".\n")
    }
  }
  suppressWarnings({
    cat(message)
    cat("\n")
    cat("Print the moving standard deviation plot: \n")
    invisible(print(p))
  })
}

find_divergence <- function(x, index) {
  
  x<- x$diagnostic_summary()$num_divergent
  positions <- which(x !=0)
  
  warnings <- list()
  
  if (length(positions) > 0) {
    warnings <- paste("The following chains had transitions ended with a divergence:", 
                      paste(x, collapse = ", "))
  }
  
  return(list(warnings = Filter(Negate(is.null), warnings), index = index))
}

find_stuck <- function(x, index, iter, window_size) {
  
  x<- x[[2]][,,1]%>%array(, dim = c(iter, 1))
  

  mov_sd <- RcppRoll::roll_sd(x, window_size, fill = NA)
  
  num_stuck_transitions <- which(mov_sd ==0)%>%length()
  # 
  warnings <- list()
  
  if (num_stuck_transitions > 0) {
    warnings <- paste("Getting stuck with", 
                      paste(num_stuck_transitions, collapse = ", "), "transitions")
  }
  
  return(list(warnings = Filter(Negate(is.null), warnings), index = index))
}



find_stuck_detail <- function(x, index, total_iter, iter_per_chain, window_size) {
  
  x_data <- x[[2]][,,1] %>% array(, dim = c(total_iter, 1))
  

  mov_sd <- RcppRoll::roll_sd(x_data, window_size, fill = NA)
  
  stuck_points <- which(mov_sd < 1.12009e-14) - (window_size / 2 - 1)
  warnings <- character(0)
  
  if (length(stuck_points) > 0) {
    
    num_chains_with_stuck <- length(unique((stuck_points - 1) %/% iter_per_chain))  # Calculate unique chains
    chains_with_stuck <- unique((stuck_points - 1) %/% iter_per_chain) + 1
    
    if (num_chains_with_stuck == 1) {
      warnings <- paste("Stuck issue is observed in 1 chain at index", chains_with_stuck)
    } else {
      warnings <- paste("Stuck issue is observed in", num_chains_with_stuck, 
                        "chains, each at chains", paste(chains_with_stuck, collapse = ", "))
      return(list(warnings = warnings, index = index, chains_with_stuck = chains_with_stuck))
    }
  }
  
}



diagnostic_ma_msd<- function(Data,  window_size, hover_threshold, hovers_lengths) {
  
  # Calculate the mean of each variable
  mean_lambda_1 <- mean(Data$lambda_1)
  mean_lambda_2 <- mean(Data$lambda_2)
  mean_lambda_3 <- mean(Data$lambda_3)
  
  # Determine which variable has the smallest mean
  if (mean_lambda_1 < mean_lambda_2 && mean_lambda_1 < mean_lambda_3) {
    smallest_mean <- "lambda_1"
  } else if (mean_lambda_2 < mean_lambda_1 && mean_lambda_2 < mean_lambda_3) {
    smallest_mean <- "lambda_2"
  } else {
    smallest_mean <- "lambda_3"
  }
  
  Data$smallest_mean<-Data[[smallest_mean]]


  mov_sd <- RcppRoll::roll_sd(Data$smallest_mean, window_size, fill = NA)
  ma <- rollmean(Data$smallest_mean, window_size,  fill = NA)
  
  plot_theme <- theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          panel.grid.major = element_line(color = "gray", linetype = "dotted"),
          panel.grid.minor = element_blank())
  
  p_msd <- ggplot(data.frame(x = seq_along(mov_sd), y = mov_sd), aes(x = x, y = y)) +
    geom_line() +
    labs(x = "Iteration", y = "Moving Standard Deviation", title = "") +
    plot_theme
  
  #p_ma <- ggplot(data.frame(x = seq_along(ma), y = ma), aes(x = x, y = y)) +
  # geom_line() +
  # labs(x = "Iteration", y = "Moving Average", title = "") +
  # plot_theme
  
  hover_points <- which(
    c(FALSE, diff(abs(mov_sd[!is.na(mov_sd)])) < hover_threshold) & 
      c(FALSE, diff(ma[!is.na(ma)])) < hover_threshold & 
      ma[!is.na(ma)] < 0.1
  )
  
  
  detect_hover_periods <- function(points) {
    starts <- numeric()
    lengths <- numeric()
    if (length(points) > 0) {
      starts <- points[1]
      current_length <- 1 + window_size -1
      for (j in 2:length(points)) {
        if (points[j] == points[j - 1] + 1) {
          current_length <- current_length + 1
        } else {
          lengths <- c(lengths, current_length)
          starts <- c(starts, points[j])
          current_length <- 1 + window_size -1
        }
      }
      lengths <- c(lengths, current_length)
    }
    list(starts = starts, lengths = lengths)
  }
  
  hovers <- detect_hover_periods(hover_points)
  
  hover_starts_long <- hovers$starts[which(hovers$lengths > hovers_lengths)] #- (window_size / 2 - 1) because we do not include NAs
  
  warnings <- if (length(hover_starts_long) > 0) {
    paste("The chain hovered around nearly zero approximately", length(hover_starts_long),
          "times, starting at iterations:", paste(hover_starts_long, collapse = ", "),
          ", with each g lasting more than ", hovers_lengths, " iterations.")
  } else {
    "No significant hovering detected."
  }
  suppressWarnings({
    cat(warnings, "\n")
    cat("Print the moving standard deviation plot: \n")

    #plot_grid <- grid.arrange(p_msd, p_ma, ncol=1)
    #grid.draw(plot_grid)
    invisible(print(p_msd))
  })
  #return(list(hovers_starts<-hovers$starts+ (window_size / 2 - 1),
  #            hovers_lengths <- hovers$lengths))
}

# Define a custom theme
my_custom_theme <- function() {
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.justification = "left"
  )
}

stuck_by_chain <- function(x, index, total_iter, iter_per_chain, window_size, stuck_length) {
  
  x_data <- x[[2]][,,1] %>% array(, dim = c(total_iter, 1))
  

  mov_sd <- RcppRoll::roll_sd(x_data, window_size, fill = NA)
  
  stuck_points <- which(mov_sd < 1.12009e-14) - (window_size / 2 - 1)
  warnings <- character(0)
  stuck_starts <- c()
  stuck_lengths <- c()
  
  if (length(stuck_points) > 0) {
    first_stuck_point <- stuck_points[1]
    stuck_starts[1] <- first_stuck_point
    current_length <- 1 + (window_size  - 1)
    
    for (j in 2:length(stuck_points)) {
      if (stuck_points[j] == stuck_points[j - 1] + 1) {
        current_length <- current_length + 1
      } else {
        stuck_lengths <- c(stuck_lengths, current_length)
        stuck_starts <- c(stuck_starts, stuck_points[j])
        current_length <- 1 + (window_size  - 1)
      }
    }
    stuck_lengths <- c(stuck_lengths, current_length)
  }
  
  stuck_points <- stuck_starts[which(stuck_lengths > stuck_length)]
  num_chains_with_stuck <- length(unique((stuck_points - 1) %/% iter_per_chain)) 
  chains_with_stuck <- unique((stuck_points - 1) %/% iter_per_chain) + 1
  persistent_stuck_chain <- (stuck_starts[which(stuck_lengths==iter_per_chain)]-1)/1000 + 1
  if (length(stuck_points) > 0) {
    
    #num_chains_with_stuck <- length(unique((stuck_points - 1) %/% iter_per_chain))  # Calculate unique chains
    #chains_with_stuck <- unique((stuck_points - 1) %/% iter_per_chain) + 1
    
    if (num_chains_with_stuck == 1) {
      warnings <- paste("Stuck issue is observed in 1 chain at index", chains_with_stuck)
    } else {
      warnings <- paste("Stuck issue is observed in", num_chains_with_stuck, 
                        "chains, each at chains", paste(chains_with_stuck, collapse = ", "))
    }
  }
  return(list(warnings = warnings, index = index, 
              num_chains_with_stuck = num_chains_with_stuck,
              chains_with_stuck = chains_with_stuck,
              persistent_stuck_chain = persistent_stuck_chain,
              stuck_lengths = stuck_lengths))
}


split_data_into_arrays <- function(data, chains_per_array) {
  # Original dimensions
  original_dimensions <- dim(data)
  # Calculate the number of resulting arrays
  num_arrays <- original_dimensions[2] / chains_per_array
  
  # Create a list to store the arrays
  array_list <- vector("list", length = num_arrays)
  
  # Loop to split the data into arrays
  for (i in 1:num_arrays) {
    start_col <- (i - 1) * chains_per_array + 1
    end_col <- i * chains_per_array
    array_list[[i]] <- data[, start_col:end_col, , drop = FALSE]
  }
  
  return(array_list)
}
Rhat_diag_by_chains <- function(resulting_arrays, small_threshold) {
  num_arrays <- length(resulting_arrays)
  Rhat_mean_by_array <- sapply(resulting_arrays, function(array) {
    rhat_values <- summarise_draws(array)$rhat
    large_rhat_positions <- which(rhat_values > small_threshold)
    larger_rhat_positions <- which(rhat_values > 1.1)
    no_large_rhat <- length(large_rhat_positions) == 0
    no_larger_rhat <- length(larger_rhat_positions) == 0
    
    return(mean(rhat_values))
  })
  
  no_large_rhat_array <- which(sapply(resulting_arrays, function(array) {
    rhat_values <- summarise_draws(array)$rhat
    length(which(rhat_values > small_threshold)) == 0
  }) == TRUE)
  
  no_larger_rhat_array <- which(sapply(resulting_arrays, function(array) {
    rhat_values <- summarise_draws(array)$rhat
    length(which(rhat_values > 1.1)) == 0
  }) == TRUE)
  
  #message("Among the ", num_arrays, " resulting chain arrays:")
  #message("- ", length(no_large_rhat_array), " arrays have Rhat values below 1.05.")
  #message("- ", length(no_larger_rhat_array), " arrays have Rhat values below 1.10.")
  
  message("Out of ", length(resulting_arrays), " 4-chain batches, ", 
          length(resulting_arrays) - length(no_larger_rhat_array), " (",
          round((length(resulting_arrays) - length(no_larger_rhat_array))/length(resulting_arrays),4)*100, "%) ",
          "have parameters with Rhat values greater than 1.10.")
  #message("The average Rhat across all parameters, by array, is:",
  #        paste(round(Rhat_mean_by_array, 3), collapse = " "))
  
  return(list(Rhat_mean_by_array = Rhat_mean_by_array))
}


fluctuating_diag_by_chain <- function(Data, window_size, selected_chains, iter_per_chain, top_percentile_threshold){
  Data1 <- Data[[1]]
  # Determine the starting and ending indices for each chain in selected_chains
  indices <- lapply(selected_chains, function(chain) {
    start_index <- (chain - 1) * iter_per_chain + 1
    end_index <- chain * iter_per_chain
    return(start_index:end_index)
  })
  
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  
  subset_data <- Data1[subset_indices, ]
  
  
  # Calculate moving standard deviation and rolling mean
  mov_sd <- RcppRoll::roll_sd(subset_data$lambda_1, window_size, fill = NA)
  
  
  ####################### Miniscule ####################
  ma <- RcppRoll::roll_mean(subset_data$lambda_1, window_size, fill = NA)
  data_matrix <- as.matrix(cbind(na.omit(ma), na.omit(mov_sd)))
  
  # Run K-means clustering
  kmeans_result <- kmeans(data_matrix, centers = 2)
  
  # Extract K-means clustering results
  centers <- kmeans_result$centers
  withinss <- kmeans_result$withinss
  
  # Create a data frame with cluster centers and withinss
  cluster_info <- data.frame(cluster = 1:nrow(centers), center_ma = centers[, 1], center_mov_sd = centers[, 2], withinss = withinss)
  
  # Find the cluster with the smallest means in 'ma' and 'mov_sd' and the smallest withinss
  smallest_cluster <- cluster_info[which.min(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  largest_cluster <- cluster_info[which.max(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  
  # Create a data frame for plotting
  data <- data.frame(`moving average` = na.omit(ma), `moving sd` = na.omit(mov_sd), 
                     cluster = kmeans_result$cluster, index = 1:nrow(data_matrix))
  # Sort the data frame by cluster
  data <- data[order(data$cluster), ]
  # Initialize a vector to store cluster-specific IDs
  cluster_ids <- numeric(0)
  
  # Generate unique IDs within each cluster
  for (i in unique(data$cluster)) {
    cluster_points <- data[data$cluster == i, ]
    num_points <- nrow(cluster_points)
    cluster_ids <- c(cluster_ids, 1:num_points)
  }
  
  
  # Add the cluster-specific IDs to the sorted data frame
  data$cluster_id <- cluster_ids
  
  k<-2
  # Calculate distances between each point and centroid
  distances <- matrix(NA, nrow = nrow(data), ncol = k)
  for (i in 1:k) {
    distances[, i] <-  (data$moving.average - kmeans_result$centers[i, 1])
  }
  
  # Add distances to the original data frame
  data2 <- cbind(data, distances)
  # Give meaningful names to distance columns
  distance_col_names <- paste0("distance_to_centroid", 1:k)
  
  # Add distances with meaningful column names to the original data frame
  data2 <- cbind(data2, setNames(distances, distance_col_names))
  colnames(data2)[(ncol(data2) - k + 1):(ncol(data2) - k + 2)] <- c("distance_to_centroid_1", "distance_to_centroid_2")
  
  
  smallest_cluster_data <- data2[data2$cluster == smallest_cluster$cluster, ]
  # Identify the top 5% largest distances for the smallest cluster
  top_5_percentile_threshold <- quantile(
    smallest_cluster_data[[paste0("distance_to_centroid_", smallest_cluster$cluster)]],
    top_percentile_threshold
  )
  
  
  # Define the distance column name
  distance_col <- paste0("distance_to_centroid_", smallest_cluster$cluster)
  
  # Add a new column indicating whether distance_to_centroid_2 is greater than top_5_percentile_threshold
  data2$above_threshold <- ifelse(
    data2$cluster == smallest_cluster$cluster &
      data2[[distance_col]] > top_5_percentile_threshold,
    "Yes",
    "No"
  )
  
  data_long <- gather(data2, key = "variable", value = "value", moving.average, moving.sd)
  
  data_long <- data_long %>%
    mutate(hover_status = ifelse(cluster == smallest_cluster$cluster, "miniscule", "not miniscule"))
  
  # Create the K-means clustering plot
  p_ma_msd <- ggplot(data_long, aes(x = index, 
                                    color = ifelse(above_threshold == "Yes", "hard to identify", hover_status))) +
    geom_point(aes(y = value), size = 0.5) +
    facet_wrap(variable ~ ., nrow = 2) +
    ggtitle(" ") +
    ylab(" ") +
    theme_minimal() + 
    labs(color = "Miniscule class") + 
    labs(x = "Iterations") +  # Set x-axis label
    my_custom_theme() +
    theme(strip.text = element_text(size = 12),
          strip.background = element_blank(),
          #axis.text.y = element_text(margin = margin(r = 50)),
          strip.placement = "outside")  +
    scale_color_manual(values = c("miniscule" = "black", 
                                  "hard to identify" = "#85c1e9", 
                                  "not miniscule" = "#FFB6C1"))
  
  
  p_trace <- ggplot(subset_data, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( expression(paste(lambda^(3))), 
                                   expression(paste(lambda^(2))), 
                                   expression(paste(lambda^(1)))))
  
  
  ###########################
  hover_points <- data$index[which(data2$cluster == smallest_cluster$cluster & 
                                     data2$above_threshold == "No")]
  sequential_positions <- diff(hover_points) == 1
  consecutive_lengths <- sequence(rle(sequential_positions)$lengths) 
  consecutive_lengths_wholeChain <- which(consecutive_lengths > 900 & sequential_positions =="TRUE")
  
  if (smallest_cluster$center_ma > 0.07 & consecutive_lengths_wholeChain%>%length() ==0) {
    warnings <- "No significant miniscule detected."
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    print(plot_grid(p_trace, p_ma_msd,  ncol = 1, align = "v"))
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  } else if (all(c(
    smallest_cluster$center_ma < 0.05,
    largest_cluster$center_ma < 0.05,
    smallest_cluster$center_mov_sd < 0.05,
    largest_cluster$center_mov_sd < 0.05
  ))) {
    warnings <- "A persistent miniscule found"
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    
    print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  }
  else
  {
    
    num_chains_with_hover <- length(unique((hover_points + 9 -1) %/% iter_per_chain))  # Calculate unique chains
    chains_with_hover <- unique((hover_points + 9 -1) %/% iter_per_chain) +1 + selected_chains[1] -1
    
    warnings <- paste("Miniscule issue is observed in", num_chains_with_hover , 
                      "chains, each at chains", paste(chains_with_hover, collapse = ", "))
    suppressWarnings({
      cat(warnings, "\n")
      cat("Print the traceplot and moving standard deviation plot: \n")
      print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
      
    })
    return(list(warnings = warnings,num_chains_with_hover = num_chains_with_hover, chains_with_hover = chains_with_hover))
  }
  
}


Rhat_diag_by_chain <- function(Data_reordered_nonPermu, num_chains, small_threshold){
  
  position <- summarise_draws(Data_reordered_nonPermu[,num_chains,])$rhat
  # Identify the positions with Rhat values greater than the threshold
  large_rhat_positions <- which(position > small_threshold)
  larger_rhat_positions <- which(position > 1.1)
  
  # Initialize an empty list to store the warning messages and positions
  warnings <- c()
  positions <- list()
  
  # If there are positions with large Rhat values
  if (length(large_rhat_positions) > 0) {
    
    # Create a warning message
    warnings$large_rhat <- paste("The following parameter positions had a Rhat greater than", small_threshold, ":", 
                                 paste(large_rhat_positions, collapse = ", "))
    
    # Store the positions in the list
    positions$large_rhat <- large_rhat_positions
    
  } else {
    # If there are no positions with large Rhat values, print a message
    print(paste("Rhat calculated out of", length(num_chains), "chains shows that NO parameters had an Rhat greater than", small_threshold, "."))
  }
  
  # Print the positions in the list
  if (length(positions$large_rhat) > 0) {
    
    message("Rhat was calculated across ", length(num_chains), " chains.  Out of 21 parameters, ", 
            length(large_rhat_positions), " had Rhat values greater than ", small_threshold, ". 
            The average Rhat across all parameters is ", round(mean(position),3), " with ", 
            length(larger_rhat_positions), " parameters exceeding an Rhat of 1.10.")
  }
}




DI <- function(l, priors, num_classes, num_subjects, selected_chains, total_chains) {
  # Extract the mcmc draws
  #mcmc_draws <- priors[[1]][[2]][, selected_chains, -c(1:30)]
  mcmc_draws <- l[, selected_chains, ]
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


twinlike_classes <- function(DI_data, selected_chains, total_chains, iter_per_chain, high_DI_value, persist_length,  happen_times) {
  Data <- DI_data
  DI_plot <- ggplot(Data, aes(x = x)) +
    geom_line(aes(y = Data[,2], color = paste(colnames(Data)[2]), linetype = paste(colnames(Data)[2])), linewidth = 1) +
    geom_line(aes(y = Data[,3], color = paste(colnames(Data)[3]), linetype = paste(colnames(Data)[3])), linewidth = 1) +
    geom_line(aes(y = Data[,4], color = paste(colnames(Data)[4]), linetype = paste(colnames(Data)[4])), linewidth = 1) +
    theme_minimal() +
    my_custom_theme() +
    labs(title = " ",
         x = "Iterations",
         y = "Distinguishability index",
         color = "Class pair") +
    scale_color_manual(values = c("#D87093", "#3498db", "black")) +
    scale_linetype_manual(values = c( "dashed","solid", "dotted"),
                          guide = guide_legend(title = "Class pair"))
  
  high_DI <- unique(c(which(Data[,2] > high_DI_value), which(Data[,3] > high_DI_value), which(Data[,4] > high_DI_value)))
  high_DI_starts <- c()
  high_DI_lengths <- c()
  if (length(high_DI) > 0) {
        first_high_DI_point <- high_DI[1]
        high_DI_starts[1] <- first_high_DI_point
        current_length <- 1 
        
          for (j in 2:length(high_DI)) {
                if (high_DI[j] == high_DI[j - 1] + 1) {
                      current_length <- current_length + 1
                  } else {
                        high_DI_lengths <- c(high_DI_lengths, current_length)
                        high_DI_starts <- c(high_DI_starts, high_DI[j])
                        current_length <- 1 
                   }
           }
       high_DI_lengths <- c(high_DI_lengths, current_length)
   
     }
  persist_high_DI <- high_DI_starts[which(high_DI_lengths > persist_length)]
  high_DI_chain <- unique((persist_high_DI  - 1) %/% iter_per_chain) + 1
  freq_table <- table((persist_high_DI  - 1) %/% iter_per_chain + 1)
  filtered_DI_values <- as.numeric(names(freq_table[freq_table >( happen_times - 1)]))
  # traceplot
  traceplot <- traceplot(Data_reordered=traceData(priors, 1, iter_per_chain*total_chains)$data, 
            num_chains = selected_chains, 
            iterations_per_chain=iter_per_chain)$traceplot.by.chain
  return_list<-list()
  return_list$DI_plot <- DI_plot
  return_list$traceplot <- traceplot
  return_list$DIplot_traceplot <- plot_grid(DI_plot, traceplot, ncol = 1, align = "v")
  return_list$filtered_DI_values <- filtered_DI_values
    return(return_list)

}

diagnostics_graphs <- function(Data, window_size, selected_chains, iter_per_chain, top_percentile_threshold,
                                      DI_data , total_chains, high_DI_value, persist_length,  happen_times){
  Data1 <- Data[[1]]
  # Determine the starting and ending indices for each chain in selected_chains
  indices <- lapply(selected_chains, function(chain) {
    start_index <- (chain - 1) * iter_per_chain + 1
    end_index <- chain * iter_per_chain
    return(start_index:end_index)
  })
  
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  
  subset_data <- Data1[subset_indices, ]
  
  
  # Calculate moving standard deviation and rolling mean
  mov_sd <- RcppRoll::roll_sd(subset_data$lambda_1, window_size, fill = NA)
  
  
  ####################### miniscule ####################
  ma <- RcppRoll::roll_mean(subset_data$lambda_1, window_size, fill = NA)
  data_matrix <- as.matrix(cbind(na.omit(ma), na.omit(mov_sd)))
  
  # Run K-means clustering
  kmeans_result <- kmeans(data_matrix, centers = 2)
  
  # Extract K-means clustering results
  centers <- kmeans_result$centers
  withinss <- kmeans_result$withinss
  
  # Create a data frame with cluster centers and withinss
  cluster_info <- data.frame(cluster = 1:nrow(centers), center_ma = centers[, 1], center_mov_sd = centers[, 2], withinss = withinss)
  
  # Find the cluster with the smallest means in 'ma' and 'mov_sd' and the smallest withinss
  smallest_cluster <- cluster_info[which.min(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  largest_cluster <- cluster_info[which.max(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  
  # Create a data frame for plotting
  data <- data.frame(`moving average` = na.omit(ma), `moving sd` = na.omit(mov_sd), 
                     cluster = kmeans_result$cluster, index = 1:nrow(data_matrix))
  # Sort the data frame by cluster
  data <- data[order(data$cluster), ]
  # Initialize a vector to store cluster-specific IDs
  cluster_ids <- numeric(0)
  
  # Generate unique IDs within each cluster
  for (i in unique(data$cluster)) {
    cluster_points <- data[data$cluster == i, ]
    num_points <- nrow(cluster_points)
    cluster_ids <- c(cluster_ids, 1:num_points)
  }
  
  
  # Add the cluster-specific IDs to the sorted data frame
  data$cluster_id <- cluster_ids
  
  k<-2
  # Calculate distances between each point and centroid
  distances <- matrix(NA, nrow = nrow(data), ncol = k)
  for (i in 1:k) {
    distances[, i] <-  (data$moving.average - kmeans_result$centers[i, 1])
  }
  
  # Add distances to the original data frame
  data2 <- cbind(data, distances)
  # Give meaningful names to distance columns
  distance_col_names <- paste0("distance_to_centroid", 1:k)
  
  # Add distances with meaningful column names to the original data frame
  data2 <- cbind(data2, setNames(distances, distance_col_names))
  colnames(data2)[(ncol(data2) - k + 1):(ncol(data2) - k + 2)] <- c("distance_to_centroid_1", "distance_to_centroid_2")
  
  
  smallest_cluster_data <- data2[data2$cluster == smallest_cluster$cluster, ]
  # Identify the top 5% largest distances for the smallest cluster
  top_5_percentile_threshold <- quantile(
    smallest_cluster_data[[paste0("distance_to_centroid_", smallest_cluster$cluster)]],
    top_percentile_threshold
  )
  
  
  # Define the distance column name
  distance_col <- paste0("distance_to_centroid_", smallest_cluster$cluster)
  
  # Add a new column indicating whether distance_to_centroid_2 is greater than top_5_percentile_threshold
  data2$above_threshold <- ifelse(
    data2$cluster == smallest_cluster$cluster &
      data2[[distance_col]] > top_5_percentile_threshold,
    "Yes",
    "No"
  )
  
  data_long <- gather(data2, key = "variable", value = "value", moving.average, moving.sd)
  
  data_long <- data_long %>%
    mutate(hover_status = ifelse(cluster == smallest_cluster$cluster, "miniscule", "not miniscule"))
  
  # Create the K-means clustering plot
  p_ma_msd <- ggplot(data_long, aes(x = index, 
                                    color = ifelse(above_threshold == "Yes", "hard to identify", hover_status))) +
    geom_point(aes(y = value), size = 0.5) +
    facet_wrap(variable ~ ., nrow = 2) +
    ggtitle(" ") +
    ylab(" ") +
    theme_minimal() + 
    labs(color = "Miniscule class") + 
    labs(x = "Iterations") +  # Set x-axis label
    my_custom_theme() +
    theme(strip.text = element_text(size = 12),
          strip.background = element_blank(),
          #axis.text.y = element_text(margin = margin(r = 50)),
          strip.placement = "outside")  +
    scale_color_manual(values = c("miniscule" = "black", 
                                  "hard to identify" = "#85c1e9", 
                                  "not miniscule" = "#FFB6C1"))
  
  
  p_trace <- ggplot(subset_data, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( "Lambda 3" = expression(paste(lambda^(3))), 
                                   "Lambda 2" = expression(paste(lambda^(2))), 
                                   "Lambda 1" = expression(paste(lambda^(1)))))
  ##############################
  source("C:/Users/xingy/OneDrive/2023 Fall/BayesIdentify/twinlike_diag.R")
  DI_plot <- twinlike_classes(DI_data, selected_chains, total_chains, iter_per_chain, 
                              high_DI_value, persist_length,  happen_times)$DI_plot
  filtered_DI_values <- twinlike_classes(DI_data, selected_chains, total_chains, iter_per_chain, 
                                         high_DI_value, persist_length,  happen_times)$filtered_DI_values
  ###########################
  hover_points <- data$index[which(data2$cluster == smallest_cluster$cluster & 
                                     data2$above_threshold == "No")]
  sequential_positions <- diff(hover_points) == 1
  consecutive_lengths <- sequence(rle(sequential_positions)$lengths) 
  consecutive_lengths_wholeChain <- which(consecutive_lengths > 900 & sequential_positions =="TRUE")
  plots <- plot_grid(p_trace, p_ma_msd,DI_plot, ncol = 1, align = "v")
  if (smallest_cluster$center_ma > 0.07 & consecutive_lengths_wholeChain%>%length() ==0) {
    warnings <- "No significant miniscule detected."
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    print(plots)
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  } else if (all(c(
    smallest_cluster$center_ma < 0.05,
    largest_cluster$center_ma < 0.05,
    smallest_cluster$center_mov_sd < 0.05,
    largest_cluster$center_mov_sd < 0.05
  ))) {
    warnings <- "A persistent miniscule found"
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    
    print(plots)
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  }
  else
  {
    
    num_chains_with_hover <- length(unique((hover_points + 9 -1) %/% iter_per_chain))  # Calculate unique chains
    chains_with_hover <- unique((hover_points + 9 -1) %/% iter_per_chain) +1 + selected_chains[1] -1
    
    warnings <- paste("miniscule issue is observed in", num_chains_with_hover , 
                      "chains, each at chains", paste(chains_with_hover, collapse = ", "))
    suppressWarnings({
      cat(warnings, "\n")
      cat("Print the traceplot and moving standard deviation plot: \n")
      print(plots)
      
    })
    return(list(warnings = warnings,num_chains_with_hover = num_chains_with_hover, chains_with_hover = chains_with_hover,
                filtered_DI_values = filtered_DI_values,
                p_trace = p_trace, p_ma_msd = p_ma_msd, DI_plot = DI_plot, plots = plots))
  }
  
}

#####################

association_by_iter <- function(Data, window_size, top_percentile_threshold, DI_3_1,DI_2_1, Stan_contents_rds,
                                num_chains, iter_per_chain){
  Data1 <- Data[[1]]

  # Calculate moving standard deviation and rolling mean
  mov_sd <- RcppRoll::roll_sd(Data1$lambda_1, window_size, fill = NA)
  
  
  ####################### miniscule ####################
  ma <- RcppRoll::roll_mean(Data1$lambda_1, window_size, fill = NA)
  data_matrix <- as.matrix(cbind(na.omit(ma), na.omit(mov_sd)))
  
  # Run K-means clustering
  kmeans_result <- kmeans(data_matrix, centers = 2)
  
  # Extract K-means clustering results
  centers <- kmeans_result$centers
  withinss <- kmeans_result$withinss
  
  # Create a data frame with cluster centers and withinss
  cluster_info <- data.frame(cluster = 1:nrow(centers), center_ma = centers[, 1], center_mov_sd = centers[, 2], withinss = withinss)
  
  # Find the cluster with the smallest means in 'ma' and 'mov_sd' and the smallest withinss
  smallest_cluster <- cluster_info[which.min(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  largest_cluster <- cluster_info[which.max(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  
  # Create a data frame for plotting
  data <- data.frame(`moving average` = na.omit(ma), `moving sd` = na.omit(mov_sd), 
                     cluster = kmeans_result$cluster, index = 1:nrow(data_matrix))
  # Sort the data frame by cluster
  data <- data[order(data$cluster), ]
  # Initialize a vector to store cluster-specific IDs
  cluster_ids <- numeric(0)
  
  # Generate unique IDs within each cluster
  for (i in unique(data$cluster)) {
    cluster_points <- data[data$cluster == i, ]
    num_points <- nrow(cluster_points)
    cluster_ids <- c(cluster_ids, 1:num_points)
  }
  
  
  # Add the cluster-specific IDs to the sorted data frame
  data$cluster_id <- cluster_ids
  
  k<-2
  # Calculate distances between each point and centroid
  distances <- matrix(NA, nrow = nrow(data), ncol = k)
  for (i in 1:k) {
    distances[, i] <-  (data$moving.average - kmeans_result$centers[i, 1])
  }
  
  # Add distances to the original data frame
  data2 <- cbind(data, distances)
  # Give meaningful names to distance columns
  distance_col_names <- paste0("distance_to_centroid", 1:k)
  
  # Add distances with meaningful column names to the original data frame
  data2 <- cbind(data2, setNames(distances, distance_col_names))
  colnames(data2)[(ncol(data2) - k + 1):(ncol(data2) - k + 2)] <- c("distance_to_centroid_1", "distance_to_centroid_2")
  
  
  smallest_cluster_data <- data2[data2$cluster == smallest_cluster$cluster, ]
  # Identify the top 5% largest distances for the smallest cluster
  top_5_percentile_threshold <- quantile(
    smallest_cluster_data[[paste0("distance_to_centroid_", smallest_cluster$cluster)]],
    top_percentile_threshold
  )
  
  
  # Define the distance column name
  distance_col <- paste0("distance_to_centroid_", smallest_cluster$cluster)
  
  # Add a new column indicating whether distance_to_centroid_2 is greater than top_5_percentile_threshold
  data2$above_threshold <- ifelse(
    data2$cluster == smallest_cluster$cluster &
      data2[[distance_col]] > top_5_percentile_threshold,
    "Yes",
    "No"
  )
  
  # Function to add missing values to the specified indices
  add_missing_values <- function(x, window_size) {
    x <- c(rep(NA, times = (window_size/2-1)), x, rep(NA, times = (window_size/2)))
    return(x)
  }
  
  # Apply the function to add missing values to the specified columns
  data_k_means <- data.frame(
    above_threshold = add_missing_values(data2$above_threshold, window_size),
    cluster = add_missing_values(data2$cluster, window_size),
    moving.average = add_missing_values(data2$moving.average, window_size),
    moving.sd = add_missing_values(data2$moving.sd, window_size)
  )
  
  data_k_means <- data_k_means %>%
    mutate(
      hover_status = ifelse(cluster == smallest_cluster$cluster, "1", "0"),
      hover_status = ifelse(above_threshold == "Yes", "0.5", hover_status) # 1.5 is hard_to_identify
    )
  dat<-nuts_params(Stan_contents_rds[[1]])%>%as.data.frame()
  accept_stat<-dat[dat$Parameter =="accept_stat__",]
  treedepth<-dat[dat$Parameter =="treedepth__",]
  divergent<-dat[dat$Parameter =="divergent__",]
  energy<-dat[dat$Parameter =="energy__",]
  stepsize<-dat[dat$Parameter =="stepsize__",]
  n_leapfrog<-dat[dat$Parameter =="n_leapfrog__",]
  
  return(data_assoc <- data.frame(chains = rep(1:num_chains, times = 1, each = iter_per_chain),
                                  iterations = rep(1:iter_per_chain, times = num_chains, each = 1),
                                  lambda_1 = Data1$lambda_1, 
                                  DI_3_1 = DI_3_1, 
                                  DI_2_1 = DI_2_1, 
                                  K_means_hover = data_k_means$hover_status,
                                  accept_stat = accept_stat$Value,
                                  treedepth = treedepth$Value,
                                  divergent = divergent$Value,
                                  energy = energy$Value,
                                  stepsize = stepsize$Value,
                                  n_leapfrog = n_leapfrog$Value))
}

hovering_class <- function(Data, window_size, iter_per_chain, top_percentile_threshold) {
  Data1 <- Data[[1]]
  
  
  # Calculate moving standard deviation and rolling mean
  mov_sd <- RcppRoll::roll_sd(Data1$lambda_1, window_size, fill = NA)
  
  ####################### miniscule ####################
  ma <- RcppRoll::roll_mean(Data1$lambda_1, window_size, fill = NA)
  data_matrix <- as.matrix(cbind(na.omit(ma), na.omit(mov_sd)))
  
  # Run K-means clustering
  kmeans_result <- kmeans(data_matrix, centers = 2)
  
  # Extract K-means clustering results
  centers <- kmeans_result$centers
  withinss <- kmeans_result$withinss
  
  # Create a data frame with cluster centers and withinss
  cluster_info <- data.frame(cluster = 1:nrow(centers), center_ma = centers[, 1], center_mov_sd = centers[, 2], withinss = withinss)
  
  # Find the cluster with the smallest means in 'ma' and 'mov_sd' and the smallest withinss
  smallest_cluster <- cluster_info[which.min(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  largest_cluster <- cluster_info[which.max(cluster_info$center_ma + cluster_info$center_mov_sd ), ]
  
  
  # Create a data frame for plotting
  data <- data.frame(`moving average` = na.omit(ma), `moving sd` = na.omit(mov_sd), cluster = kmeans_result$cluster, index = 1:nrow(data_matrix))
  # Sort the data frame by cluster
  data <- data[order(data$cluster), ]
  # Initialize a vector to store cluster-specific IDs
  cluster_ids <- numeric(0)
  
  # Generate unique IDs within each cluster
  for (i in unique(data$cluster)) {
    cluster_points <- data[data$cluster == i, ]
    num_points <- nrow(cluster_points)
    cluster_ids <- c(cluster_ids, 1:num_points)
  }
  
  # Add the cluster-specific IDs to the sorted data frame
  data$cluster_id <- cluster_ids
  
  
  
  k<-2
  # Calculate distances between each point and centroid
  distances <- matrix(NA, nrow = nrow(data), ncol = k)
  for (i in 1:k) {
    distances[, i] <-  (data$moving.average - kmeans_result$centers[i, 1])
  }
  
  
  # Add distances to the original data frame
  data2 <- cbind(data, distances)
  # Give meaningful names to distance columns
  distance_col_names <- paste0("distance_to_centroid", 1:k)
  
  # Add distances with meaningful column names to the original data frame
  data2 <- cbind(data2, setNames(distances, distance_col_names))
  colnames(data2)[(ncol(data2) - k + 1):(ncol(data2) - k + 2)] <- c("distance_to_centroid_1", "distance_to_centroid_2")
  
  
  smallest_cluster_data <- data2[data2$cluster == smallest_cluster$cluster, ]
  # Identify the top 5% largest distances for the smallest cluster
  top_5_percentile_threshold <- quantile(
    smallest_cluster_data[[paste0("distance_to_centroid_", smallest_cluster$cluster)]],
    top_percentile_threshold
  )
  
  
  # Define the distance column name
  distance_col <- paste0("distance_to_centroid_", smallest_cluster$cluster)
  
  # Add a new column indicating whether distance_to_centroid_2 is greater than top_5_percentile_threshold
  data2$above_threshold <- ifelse(
    data2$cluster == smallest_cluster$cluster &
      data2[[distance_col]] > top_5_percentile_threshold,
    "Yes",
    "No"
  )
  
  data_long <- gather(data2, key = "variable", value = "value", moving.average, moving.sd)
  
  data_long <- data_long %>%
    mutate(hover_status = ifelse(cluster == smallest_cluster$cluster, "miniscule", "not miniscule"))
  
  # Create the K-means clustering plot
  p_ma_msd <- ggplot(data_long, aes(x = index, 
                                    color = ifelse(above_threshold == "Yes", "hard to identify", hover_status))) +
    geom_point(aes(y = value), size = 0.5) +
    facet_wrap(variable ~ ., nrow = 2) +
    ggtitle(" ") +
    ylab(" ") +
    theme_minimal() + 
    labs(color = "Miniscule class") + 
    labs(x = "Iterations") +  # Set x-axis label
    my_custom_theme() +
    theme(strip.text = element_text(size = 12),
          strip.background = element_blank(),
          #axis.text.y = element_text(margin = margin(r = 50)),
          strip.placement = "outside")  +
    scale_color_manual(values = c("miniscule" = "black", 
                                  "hard to identify" = "#85c1e9", 
                                  "not miniscule" = "#FFB6C1"))
  
  
  # traceplot
  p_trace <- Data$traceplot
  
  # two distinguishable classes
  diff32 <-abs((Data1$mu_intercept_3 + Data1$mu_slope_3 * 1.5 + Data1$sq_mu_slope_3 * 2.25)-
                 (Data1$mu_intercept_2 + Data1$mu_slope_2 * 1.5 + Data1$sq_mu_slope_2 * 2.25))
  diff31 <- abs((Data1$mu_intercept_3 + Data1$mu_slope_3 * 1.5 + Data1$sq_mu_slope_3 * 2.25)-
                  (Data1$mu_intercept_1 + Data1$mu_slope_1 * 1.5 + Data1$sq_mu_slope_1 * 2.25))
  diff12 <- abs((Data1$mu_intercept_1 + Data1$mu_slope_1 * 1.5 + Data1$sq_mu_slope_1 * 2.25)-
                  (Data1$mu_intercept_2 + Data1$mu_slope_2 * 1.5 + Data1$sq_mu_slope_2 * 2.25))
  
  xseq <- seq(0, 3, length.out = 100)
  
  
  # Create quadratic regression lines
  yseq3 <- Data1$mu_intercept_3 + Data1$mu_slope_3 * xseq + Data1$sq_mu_slope_3 * xseq^2
  yseq1 <- Data1$mu_intercept_1 + Data1$mu_slope_1 * xseq + Data1$sq_mu_slope_1 * xseq^2
  yseq2 <- Data1$mu_intercept_2 + Data1$mu_slope_2 * xseq + Data1$sq_mu_slope_2 * xseq^2
  
  ma_diff32 <- RcppRoll::roll_mean(diff32, window_size, fill = NA)
  ma_diff31 <- RcppRoll::roll_mean(diff31, window_size, fill = NA)
  ma_diff12 <- RcppRoll::roll_mean(diff12, window_size, fill = NA)
  
  ma_mean_diff_plot <- ggplot(Data1, aes(x = x)) +
    geom_line(aes(y = ma_diff31, color = "class 3 vs 1"), linetype = "solid", linewidth = 1) +
    geom_line(aes(y = ma_diff32, color = "class 2 vs 3"), linetype = "solid", linewidth = 1) +
    geom_line(aes(y = ma_diff12, color = "class 1 vs 2"), linetype = "solid", linewidth = 1) +
    geom_hline(aes(yintercept=0),linetype = "dotted") + 
    labs(x = "Iterations") +
    # Add titles and labels
    labs(title = " ",
         x = "",
         y = "MA of mean differences when t=1.5",
         color = "Difference between classes") +
    # Adjust theme
    theme_minimal() +
    my_custom_theme()   +
    scale_color_manual(values = c("class 3 vs 1" = "#D87093" ,
                                  "class 2 vs 3" = "#3498db" , 
                                  "class 1 vs 2" = "black")) 
  
  
  if (smallest_cluster$center_ma > 0.05) {
    warnings <- "No significant miniscule detected."
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    print(plot_grid(p_trace, p_ma_msd,ma_mean_diff_plot, ncol = 1, align = "v"))
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  } else if (all(c(
    smallest_cluster$center_ma < 0.05,
    largest_cluster$center_ma < 0.05,
    smallest_cluster$center_mov_sd < 0.05,
    largest_cluster$center_mov_sd < 0.05
  ))) {
    warnings <- "A persistent miniscule found"
    cat(warnings, "\n")
    cat("Print the traceplot: \n")
    
    print(plot_grid(p_trace, p_ma_msd,ma_mean_diff_plot, ncol = 1, align = "v"))
    #print(plot_grid(p_trace, p_ma_msd, ncol = 1, align = "v"))
  }
  else
  {hover_points <- data$index[which(data2$cluster == smallest_cluster$cluster & 
                                      data2$above_threshold == "No")]
  
  
  num_chains_with_hover <- length(unique((hover_points + 9 -1) %/% iter_per_chain))  # Calculate unique chains
  chains_with_hover <- unique((hover_points + 9 -1) %/% iter_per_chain) +1
  
  warnings <- paste("miniscule issue is observed in", num_chains_with_hover , 
                    "chains, each at chains", paste(chains_with_hover, collapse = ", "))
  suppressWarnings({
    cat(warnings, "\n")
    cat("Print the traceplot and moving standard deviation plot: \n")
    print(plot_grid(p_trace, p_ma_msd,ma_mean_diff_plot, ncol = 1, align = "v"))
    
  })
  return(list(warnings = warnings,num_chains_with_hover = num_chains_with_hover, chains_with_hover = chains_with_hover))
  }
  
  
}


traceplot <- function(Data_reordered, num_chains, iterations_per_chain){
  
  return_list<-list()
  return_list$traceplot<- ggplot(Data_reordered, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( "Lambda 3" = expression(paste(lambda^(3))), 
                                   "Lambda 2" = expression(paste(lambda^(2))), 
                                   "Lambda 1" = expression(paste(lambda^(1)))))
  
  
  # Determine the starting and ending indices for each chain in num_chains
  indices <- lapply(num_chains, function(chain) {
    start_index <- (chain - 1) * iterations_per_chain + 1
    end_index <- chain * iterations_per_chain
    return(start_index:end_index)
  })
  
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  subset_data <- Data_reordered[subset_indices, ]
  
  return_list$traceplot.by.chain <- ggplot(subset_data, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( "Lambda 3" = expression(paste(lambda^(3))), 
                                   "Lambda 2" = expression(paste(lambda^(2))), 
                                   "Lambda 1" = expression(paste(lambda^(1)))))
  
  
  return(return_list)
}

traceplot_label <- function(Data_reordered, iterations_to_display, label){
  
  return_list<-list()
  return_list$traceplot<- ggplot(Data_reordered, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    geom_vline(xintercept = label, linetype = "dashed", color = "red", linewidth = 1) +  # Add vertical line
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( "Lambda 3" = expression(paste(lambda^(3))), 
                                   "Lambda 2" = expression(paste(lambda^(2))), 
                                   "Lambda 1" = expression(paste(lambda^(1)))))
  
  
  
  subset_data <- Data_reordered[iterations_to_display, ]
  
  return_list$traceplot.by.chain <- ggplot(subset_data, aes(x=x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", linewidth = 1) +
    geom_vline(xintercept = label, linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal() +
    my_custom_theme()+
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c( "Lambda 3" = "darkred", 
                                   "Lambda 2" = "steelblue", 
                                   "Lambda 1" = "black"),
                       labels = c( "Lambda 3" = expression(paste(lambda^(3))), 
                                   "Lambda 2" = expression(paste(lambda^(2))), 
                                   "Lambda 1" = expression(paste(lambda^(1)))))
  
  
  return(return_list)
}


traceplot_stan_diag <- function(Data_reordered, Stan_contents_rds, num_chains, iterations_per_chain, stan_diag, diag_value){
  dat<-nuts_params(Stan_contents_rds[[1]])%>%as.data.frame()
  accept_stat<-dat[dat$Parameter =="accept_stat__",]
  treedepth<-dat[dat$Parameter =="treedepth__",]
  divergent<-dat[dat$Parameter =="divergent__",]
  energy<-dat[dat$Parameter =="energy__",]
  stepsize<-dat[dat$Parameter =="stepsize__",]
  n_leapfrog<-dat[dat$Parameter =="n_leapfrog__",]
  return_list<-list()
  
  
  
  # Determine the starting and ending indices for each chain in num_chains
  indices <- lapply(num_chains, function(chain) {
    start_index <- (chain - 1) * iterations_per_chain + 1
    end_index <- chain * iterations_per_chain
    return(start_index:end_index)
  })
  
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  subset_data <- Data_reordered[subset_indices, ]
  
  return_list$traceplot.by.chain <- ggplot(subset_data, aes(x = x)) +
    geom_line(aes(y = lambda_3, color = "Lambda 3"), linetype = "dotted", size = 1) +
    geom_line(aes(y = lambda_2, color = "Lambda 2"), linetype = "dashed", size = 1) +
    geom_line(aes(y = lambda_1, color = "Lambda 1"), linetype = "solid", size = 1) +
    geom_point(data = subset_data[which(get(paste("stan_diag"))$Value[get(paste("subset_indices"))] == diag_value), ],
               aes(x = x, y = -0.0001), color = "red", size = 0.5) +
    theme_minimal() +
    my_custom_theme() +
    labs(title = " ",
         x = "Iterations",
         y = "Class probabilities",
         color = "Class probabilities") +
    scale_color_manual(values = c("Lambda 3" = "darkred",
                                  "Lambda 2" = "steelblue",
                                  "Lambda 1" = "black"),
                       labels = c("Lambda 3" = expression(paste(lambda^(3))),
                                  "Lambda 2" = expression(paste(lambda^(2))),
                                  "Lambda 1" = expression(paste(lambda^(1)))))
  
  
  
  return(return_list)
}
