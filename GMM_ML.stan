functions {
  // Custom function for multivariate normal likelihood
  real mmn(real beta_1, real beta_2, real beta_3,
          vector time_seg, vector time_sq_seg,
          matrix Sigma, real sigma_e, vector y_seg) {
    vector[rows(time_seg)] mu_seg;
    matrix[rows(time_seg), 2] Z_j;
    matrix[rows(time_seg), rows(time_seg)] Cov_j;
    mu_seg = beta_1 + beta_2 * time_seg + beta_3 * time_sq_seg;  // Equation (1)
    Z_j = append_col(rep_vector(1.0, rows(time_seg)), time_seg);
    Cov_j = Z_j * Sigma * Z_j' + diag_matrix(rep_vector(sigma_e^2, rows(time_seg))); // Equation (8)
    return multi_normal_lpdf(y_seg | mu_seg, Cov_j);         // return multivariate normal log-likelihood
  }
}

// Specifies the data structures and variables that will be input into the model.
data {
  int<lower=1> N;                  // number of observations
  int<lower=1> J;                  // number of subjects
  array[N] int<lower=1, upper=J> Subject;   // subject ID for each observation
  array[J] int<lower=1> s;         // subject size
  vector[N] y;                     // response variable
  vector[N] time;                  // the time variable
  vector[N] time_sq;               // squared time variable
  int<lower=1> K;                  // number of latent classes
  int<lower=0> Dir_alpha;          // concentration parameter of the Dirichlet prior
  int<lower=0> Cauchy_scale;       // scale parameter of the half-Cauchy prior
  int<lower=0> Normal_scale;       // scale parameter of the half-Normal prior
  int<lower=0> prior;              // indicator for the type of prior used for standard deviation parameters
  real<lower=0> eta;                // shape parameter for the LKJ distribution
}

// Allows for transformations of the input data before it is used in the modeling process.
transformed data {
  array[J] int pos;               // position indices indicating subject boundaries
  int new_j = 1;                  // counter for new subject index
  int curr_sub = Subject[1];      // current subject identifier
  pos[1] = 1;                     // start position for the first subject
  // Determine positions where subjects change
  for (n in 1:N) {
    if (Subject[n] != curr_sub) { // check if current subject differs from previous
      new_j += 1;                 // increment new subject index
      pos[new_j] = n;             // record position where new subject starts
      curr_sub = Subject[n];      // update current subject identifier
    }
  }
}

// Declares the model parameters that Stan will estimate.
parameters {
  simplex[K] lambda;    // class probability parameter vector
  vector[K] beta_1;     // fixed intercept vector 
  vector[K] beta_2;     // fixed slope of time variable vector 
  vector[K] beta_3;     // fixed slope of squared time variable vector 
  array[K] vector<lower=0>[2] sigma_u;   //  random-intercept and random-slope sds
  array[K] cholesky_factor_corr[2] L_Omega;  // Cholesky factor of the correlation matrix
  real<lower=0> sigma_e;  // residual sd
}

// Defines parameters that are deterministic transformations of the primary parameters.
transformed parameters {
  array[K] corr_matrix[2] Omega;  // correlation matrices for random intercepts and slopes
  array[K] cov_matrix[2] Sigma;   // covariance matrices for random intercepts and slopes
  array[J] vector[K] lps;   // log posterior probabilities for each subject
  // Compute correlation and covariance matrices for each class
  for (k in 1:K) {
    Omega[k] = multiply_lower_tri_self_transpose(L_Omega[k]);   // calculate correlation matrix
    Sigma[k] = quad_form_diag(Omega[k], sigma_u[k]);           // calculate covariance matrix
  }
  // Compute log posterior probabilities for each subject
  for (j in 1:Subject[N]) {
    lps[j] = log(lambda); // initialize with log class probabilities
    for (k in 1:K) {
      lps[j][k] += mmn(beta_1[k], beta_2[k], beta_3[k],
        segment(time, pos[j], s[j]), segment(time_sq, pos[j], s[j]),
        Sigma[k], sigma_e, segment(y, pos[j], s[j])); // add log-likelihood contribution
    }
  }
}

// Specifies the likelihood function, priors, and the increment of the log probability density.
model {
  // Priors
  for (k in 1:K) {
    beta_1[k] ~ normal(0, 10);
    beta_2[k] ~ normal(0, 5);
    beta_3[k] ~ normal(0, 5);
    if (prior == 0) {
      sigma_u[k] ~ cauchy(0, Cauchy_scale);
    }
    if (prior == 1) {
      sigma_u[k] ~ normal(0, Normal_scale);
    }
    L_Omega[k] ~ lkj_corr_cholesky(eta);  
  }
  sigma_e ~ exponential(1);
  lambda ~ dirichlet(rep_vector(Dir_alpha, K));
  // Likelihood shown in Equation (9)
  for (j in 1:Subject[N]) target += log_sum_exp(lps[j]);  // accumulate log-likelihood contributions 
}

// Used to generate predictions or derived quantities after the model has been fitted.
generated quantities {
  vector[J] log_lik;                     // log-likelihood for each subject
  array[J] int<lower=1> pred_class_dis;  // predicted class assignment for subject j in class k
  array[J] simplex[K] pred_class;        // post. prob. of subject j in class k
  // Compute log-likelihood, predicted class probabilities, and assignments
  for (j in 1:Subject[N]) {
    log_lik[j] = log_sum_exp(lps[j]);   // log-likelihood contribution for subject j
    pred_class[j] = softmax(lps[j]);    // posterior class probabilities for subject j
    pred_class_dis[j] = categorical_rng(pred_class[j]);  // predicted class assignment for subject j
  }
}