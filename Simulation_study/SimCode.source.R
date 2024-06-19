
UncGMM_data <- function(n_pers, n_time, beta_int, beta_slo_time, beta_slo_time_sq,
                        sd_i, sd_s, cor_is,
                        sd_r,
                        K, lambda_K) {
  latent_group <- sample(1:K,
                         prob = lambda_K,
                         size = n_pers,
                         replace = TRUE)
  eff = map(1:n_pers, ~ mvrnorm(
    1,
    mu = c(beta_int[latent_group[.x]],
           beta_slo_time[latent_group[.x]],
           beta_slo_time_sq[latent_group[.x]]),
    Sigma = rbind(
      c((sd_i[latent_group[.x]])^2,
        cor_is[latent_group[.x]] * sd_i[latent_group[.x]] *
          sd_s[latent_group[.x]],
        0),
      c(cor_is[latent_group[.x]] * sd_i[latent_group[.x]] *
          sd_s[latent_group[.x]], (sd_s[latent_group[.x]])^2,
        0), 
      c(0,0,0)
    )
  )) %>% do.call(rbind, .)
  colnames(eff) = c("intercept", "slope_time", "slope_time_sq")
  
  
  dat = data.frame(
    ID = rep(1:n_pers, each = n_time),
    time = rep(1:n_time, times = n_pers),
    time_sq = rep(c(0,1,4,9), times = n_pers),
    latent_group = rep(latent_group, each = n_time),
    int = rep(eff[, 1], each = n_time),
    slo = rep(eff[, 2], each = n_time),
    slo_sq = rep(eff[, 3], each = n_time),
    y = NA
  )
  dat$time <- dat$time - 1
  #dat$interact <- dat$time*dat$covar
  y = with(dat,
           (int) + (slo) * time  + (slo_sq) * time_sq + 
             rnorm( n = n_pers * n_time,
                    mean = 0,
                    sd = sqrt(sd_r)
             ))
  dat$y <- y
  return(dat)
}