source("Diagnostics/Diagnostics.source.R")
# Define a function to detect stuck sequences by chain
stuck_by_chain <- function(x, index, total_iter, iter_per_chain, window_size, stuck_length) {
  
  # Extract the relevant data for analysis
  x_data <- x[[2]][,,1] %>% array(, dim = c(total_iter, 1))
  
  # Compute the moving standard deviation
  mov_sd <- RcppRoll::roll_sd(x_data, window_size, fill = NA)
  
  # Identify points where standard deviation is below a threshold
  stuck_points <- which(mov_sd < 1.12009e-14) - (window_size / 2 - 1)
  
  # Initialize variables to store results
  warnings <- character(0)
  stuck_starts <- c()
  stuck_lengths <- c()
  
  # Check for stuck sequences
  if (length(stuck_points) > 0) {
    first_stuck_point <- stuck_points[1]
    stuck_starts[1] <- first_stuck_point
    current_length <- 1 + (window_size - 1)
    
    for (j in 2:length(stuck_points)) {
      if (stuck_points[j] == stuck_points[j - 1] + 1) {
        current_length <- current_length + 1
      } else {
        stuck_lengths <- c(stuck_lengths, current_length)
        stuck_starts <- c(stuck_starts, stuck_points[j])
        current_length <- 1 + (window_size - 1)
      }
    }
    stuck_lengths <- c(stuck_lengths, current_length)
  }
  
  # Filter out sequences longer than specified stuck length
  stuck_points <- stuck_starts[which(stuck_lengths > stuck_length)]
  
  # Count number of chains with stuck sequences
  num_chains_with_stuck <- length(unique((stuck_points - 1) %/% iter_per_chain)) 
  chains_with_stuck <- unique((stuck_points - 1) %/% iter_per_chain) + 1
  
  # Calculate persistent stuck chains
  persistent_stuck_chain <- (stuck_starts[which(stuck_lengths == iter_per_chain)] - 1) / 1000 + 1
  
  # Generate warning messages based on findings
  if (length(stuck_points) > 0) {
    if (num_chains_with_stuck == 1) {
      warnings <- paste("Stuck issue observed in 1 chain at index", chains_with_stuck)
    } else {
      warnings <- paste("Stuck issue observed in", num_chains_with_stuck, 
                        "chains, each at chains", paste(chains_with_stuck, collapse = ", "))
    }
  }
  
  # Return a list of results
  return(list(warnings = warnings, index = index, 
              num_chains_with_stuck = num_chains_with_stuck,
              chains_with_stuck = chains_with_stuck,
              persistent_stuck_chain = persistent_stuck_chain,
              stuck_lengths = stuck_lengths))
}