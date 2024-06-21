source("Diagnostics/Diagnostics.source.R")
# Function to detect twinlike-class behavior using DI plot and traceplot
twinlike_classes <- function(DI_data, selected_chains, total_chains, iter_per_chain, high_DI_value, persist_length, happen_times) {
  
  # Extract DI data
  indices <- lapply(selected_chains, function(chain) {
    start_index <- (chain - 1) * iter_per_chain + 1
    end_index <- chain * iter_per_chain
    return(start_index:end_index)
  })
  
  # Combine the indices and subset the data
  subset_indices <- unlist(indices)
  Data <- DI_data[subset_indices, ]
  
  # Generate DI plot
  DI_plot <- ggplot(Data, aes(x = x)) +
    geom_line(aes(y = Data[,2], color = paste(colnames(Data)[2]), linetype = paste(colnames(Data)[2])), linewidth = 1) +
    geom_line(aes(y = Data[,3], color = paste(colnames(Data)[3]), linetype = paste(colnames(Data)[3])), linewidth = 1) +
    geom_line(aes(y = Data[,4], color = paste(colnames(Data)[4]), linetype = paste(colnames(Data)[4])), linewidth = 1) +
    theme_minimal() +                # Apply minimal theme
    my_custom_theme() +              # Apply custom theme (assumed to be defined elsewhere)
    labs(title = " ",                 # Title for the plot (can be customized)
         x = "Iterations",            # X-axis label
         y = "Distinguishability index",  # Y-axis label
         color = "Class pair") +      # Color legend label
    scale_color_manual(values = c("#D87093", "#3498db", "black")) +  # Custom color scheme
    scale_linetype_manual(values = c("dashed", "solid", "dotted"),   # Custom linetype for legend
                          guide = guide_legend(title = "Class pair"))  # Legend title
  
  # Identify iterations with high DI values
  high_DI <- unique(c(which(Data[,2] > high_DI_value), which(Data[,3] > high_DI_value), which(Data[,4] > high_DI_value)))
  high_DI_starts <- c()
  high_DI_lengths <- c()
  
  # Calculate lengths of consecutive high DI segments
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
  
  # Identify chains with persistently high DI values
  persist_high_DI <- high_DI_starts[which(high_DI_lengths > persist_length)]
  high_DI_chain <- unique((persist_high_DI - 1) %/% iter_per_chain) + 1
  
  # Count occurrences of persistently high DI values
  freq_table <- table((persist_high_DI - 1) %/% iter_per_chain + 1)
  filtered_DI_values <- as.numeric(names(freq_table[freq_table > (happen_times - 1)]))
  
  # Generate traceplot to visualize chain behavior
  traceplot <- traceplot(Data_reordered = traceData(priors, 1, iter_per_chain * total_chains)$data, 
                         num_chains = selected_chains, 
                         iterations_per_chain = iter_per_chain)$traceplot.by.chain
  
  # Return a list of outputs
  return_list <- list()
  return_list$DI_plot <- DI_plot               # Distinguishability index plot
  return_list$traceplot <- traceplot           # Traceplot of selected chains
  return_list$DIplot_traceplot <- plot_grid(DI_plot, traceplot, ncol = 1, align = "v")  # Combined plot
  return_list$filtered_DI_values <- filtered_DI_values  # DI values exceeding occurrence threshold
  
  return(return_list)
}
