# Bayesian Identification and Estimation of Growth Mixture Models
This repository provides Stan code for fitting Growth Mixture Models (GMMs) and includes diagnostics and graphs described in our paper (link will be provided upon publication).

## Abstract

The development of reading skills in children aged 6 to 14 is examined by applying Bayesian growth mixture models (GMMs) to reading recognition data from the National Longitudinal Study of Youth (NLSY). Three types (or classes) of learners are identified who differ in both their mean growth trajectories and in the nature of the within-class heterogeneity in growth trajectories. For vague or diffuse priors, Markov Chain Monte Carlo (MCMC) methods can encounter problems for GMMs due to identifiability issues that we discuss in detail. 
We review the literature on Bayesian identification and propose a viable definition for latent variable (or hierarchical Bayesian) models, such as GMMs, based on the marginal likelihood. We also give a brief introduction to Hamiltonian Monte Carlo estimation methods implemented in Stan, demonstrate that problematic behavior occurs for the NLSY data due to identifiability issues, and 
suggest diagnostics \newtext{designed specifically} for detecting such behavior. A simulation study shows that problematic behavior is quite common in GMMs, and use of diagnostics to detect the problem is therefore important. Simulations also suggest that choosing more informative priors than the vague defaults can mitigate the problems. 

## Contents

- Stan Code
- Diagnostics
- Graphs
