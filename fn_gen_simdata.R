#############################################################################
###                                                                       ###
###   Objective: Synthetic data generation for comparative analysis       ###
###                                                                       ###
###   Author: Sangyoon Yi (sayi@okstate.edu)                              ###
###           Myeongjong Kang (mkangstat@gmail.com)                       ###
###                                                                       ###
###   Key functions: gen_pred, gen_latent, gen_resp, sim_model            ###
###                                                                       ###
###   Reference: Estimation of treatment effect in clinical trials of     ###
###                         continuous endpoints with retrieved dropouts  ###
###                                                                       ###
#############################################################################

#############################################################################
###   1. function to simulate each baseline response and predictor        ###
#############################################################################

# The function below generates Y_{i,0} at baseline and binary variable for trt/placebo

# Inputs
# n: sample size
# ntrt: the total number of subjects having a treatment
# mu: mean of Normal dist'n for baseline response
# s: standard deviation of Normal dist'n for baseline response

# Output
# pmat: n by 2 matrix of simulated data
# where the first column shows the vector of initial responses at baseline
# the second column consists of binary variable that indicates treatment assignment

gen_pred <- function(n, ntrt, mu, s){
  pmat <- matrix(NA, nrow = n, ncol = 2)
  pmat[,1] <- rnorm(n, mu, s)
  pmat[,2] <- c(rep(1, ntrt), rep(0, n - ntrt))
  return(pmat)
}

#############################################################################
###   2. function to compute R_{i} for subject categorization             ###
#############################################################################

# The function below simulates R_{i} 

# Inputs
# gvec: 3-dimensional vector of gamma (intercept and the other two) for the latent variable model
# pi: the conditional probability of Q = 1 given D = 0
# pmat: n by 2 matrix simulated from the gen_pred function above 
# sd_u: the standard deviation of error in the latent variable model
# opt_std: a logical to indicate whether the baseline response to be scaled or not

# Output
# Rvec: n-dimensional vector of R_{i} 

gen_latent <- function(gvec, pi, pmat, sd_u, opt_std = T){
  n <- nrow(pmat)
  if(opt_std){
    pmat[,1] <- scale(pmat[,1])
  }
  lvec <- cbind(1, pmat)%*%gvec + rnorm(n, mean = 0, sd = sd_u)   # Compute the latent variable
  Dvec <- ifelse(lvec < 0, 1, 0)   # D_{i} = 1 if the latent variable is "smaller" than 0
  Qvec <- rep(1, n) # Qvec: a vector of Q_{i} to further distinguish RD against Missing
  Qvec[Dvec==0] <- sample(c(1,0), sum(Dvec==0), replace = T, prob = c(pi, 1 - pi))  # Q_{i} follows a Bernoulli with pi conditional on D_{i}=0  
  Rvec <- rep(1, n) # R_{i}=1 (Complete)
  Rvec[Dvec==0 & Qvec==1] <- 0 # R_{i}=0 (RD) if D_{i}=0 and Q_{i}=1
  Rvec[Dvec==0 & Qvec==0] <- -1 # R_{i}=-1 (Missing) if D_{i}=0 and Q_{i}=0
  return(Rvec)
}

#############################################################################
###   3. function to compute the observed response for completer and RD   ###
#############################################################################

# Inputs
# Rvec: n-dimensional vector of R_{i} 
# bvec: 2-dimensional vector of beta for the linear model
# delta: effect associated with either RD or Missing
# pmat: n by 2 matrix of predictors
# sd_eps: standard deviation of error in the model for latent variable

# Outputs
# obs_Y: n-dimensional vector of observed responses (zero if R_{i} = -1)
# Y: n-dimensional vector of full response (including R_{i} = -1)

gen_resp <- function(Rvec, bvec, delta, pmat, sd_eps){
  Y <- pmat%*%bvec + delta*ifelse(Rvec<=0, 1, 0) + rnorm(length(Rvec), sd = sd_eps)
  obs_Y <- rep(0, length(Y)) ; obs_Y[Rvec>=0] <- Y[Rvec>=0]
  return(list(Y = Y, obs_Y = obs_Y))
}

#############################################################################
###   4. function to integrate the above three functions for simulation   ###
#############################################################################

# Inputs
# Same as the descriptions above

# Outputs
# pmat: n by 2 matrix of predictors
# Rvec: n-dimensional vector of R_{i}
# obs_Y: n-dimensional vector of observed responses
# Y: n-dimensional vector of full response (including R_{i} = -1)

sim_model <- function(n, ntrt, pi, delta, bvec, gvec, mu, s, sd_eps, sd_u, opt_std = T){
  pmat <- gen_pred(n, ntrt, mu, s) 
  Rvec <- gen_latent(gvec, pi, pmat, sd_u, opt_std)
  res_list <- gen_resp(Rvec, bvec, delta, pmat, sd_eps)
  return(list(pmat = pmat, Rvec = Rvec, obs_Y = res_list$obs_Y, Y = res_list$Y))
}
