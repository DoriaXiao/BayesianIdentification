source("Diagnostics/Diagnostics.source.R")
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