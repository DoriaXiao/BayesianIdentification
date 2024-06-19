library(label.switching)
library(rstan)
post_processing <- function(chains, iterations, K, J, post_class, mcmc, post_class_p, post_par){

  # set of selected relabeling algorithm
  set <-
    c("STEPHENS", "PRA", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2")
  mapindex = which.max(post_par$lp__)
  ls_gmm_long <-
    label.switching(
      method = set,
      zpivot = post_class[mapindex,],
      z = post_class,
      K = K,
      prapivot = mcmc[mapindex, ,],
      constraint = 1,
      mcmc = mcmc,
      p = post_class_p
    )
  # similarity of the classification
  similarity <- ls_gmm_long$similarity
  # permuted posterior based on ECR method
  mcmc_permuted <- permute.mcmc(mcmc, ls_gmm_long$permutations$STEPHENS)
  dim(mcmc_permuted$output)
  # change dimension for each parameter defined as in the Stan code
  mcmc_permuted <-
    array(
      data = mcmc_permuted$output,
      dim = c(iterations, chains, K*J) # (iterations, chains, classes*parameters)
    )
  fit_permuted <-
    monitor(mcmc_permuted, warmup = 0,  digits_summary = 3) %>% as.data.frame()
  
  PostProcessing <- list(similarity, mcmc_permuted, fit_permuted)
  return(PostProcessing)
}