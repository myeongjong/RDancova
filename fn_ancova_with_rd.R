#############################################################################
###                                                                       ###
###   Objective:                                                          ###
###                                                                       ###
###   Author: Sangyoon Yi (sayi@okstate.edu)                              ###
###           Myeongjong Kang (mkangstat@gmail.com)                       ###
###                                                                       ###
###   Key functions:                                                      ###
###                                                                       ###
###   Reference: Estimation of treatment effect in clinical trials of     ###
###                         continuous endpoints with retrieved dropouts  ###
###                                                                       ###
#############################################################################

#############################################################################
###   1. compute the function to compute the MLE                          ###
#############################################################################

# Inputs
# pmat: n by 2 matrix of predictors
# obs_y: n dimensional vector of observed response
# Rvec: n dimensional vector of R_{i}
# max_iter: maximum number of iterations to be passed to glm function
# opt_std: a logical to indicate whether the baseline response to be scaled or not
# opt_offset: a logical to indicate whether we set the gamma_{1} (associated with the baseline) to be 1 or not

# Outputs
# pi_hat: the estimates of pi 
# res_lmod: the resulting output from linear model
# sig_hat: the estimates of \sigma_{\epsilon} 
# res_pmod: the resulting output from probit model

#pmat = dat$pmat ; obs_y = dat$obs_Y ; Rvec = dat$Rvec ; max_iter = 500 ;opt_std = T ;opt_offset = T

find_mle <- function(pmat, obs_y, Rvec, max_iter = 500, opt_std = T, opt_offset = T){
  
  # 1. Compute \hat{\pi}
  pi_hat <- sum(Rvec==0)/sum(Rvec<1)
  
  # 2. Compute \hat{\beta}, \hat{\delta} and \hat{\sigma}_{\epsilon}
  # the indices with either R_{i}=1 or R_{i}=0
  ind_not_miss <- which(Rvec>=0)
  
  # the difference between responses would be used as new response
  new_resp <- obs_y[ind_not_miss] - pmat[ind_not_miss, 1]  
  
  # set a new design matrix
  new_dmat <- cbind(ifelse(Rvec[ind_not_miss]==0, 1, 0), pmat[ind_not_miss,])
  
  lmod <- lm(new_resp ~ new_dmat)
  sum_lmod <- summary(lmod) ; sig_hat <- sigma(lmod)
  res_lmod <- coef(sum_lmod)
  rownames(res_lmod) <- c("intercept", "delta", paste("beta", sep = " ", c(1:ncol(pmat))))
  resid_vec <- as.vector(residuals(lmod))
  
  # 3. Compute \hat{\gamma}
  Dvec <- ifelse(Rvec==1, 1, 0)  
  Dvec2 <- ifelse(Dvec==1, 0, 1) # Should we change the response since we changed the inequality for the latent variable?
  if(opt_std & opt_offset){
    pmod <- glm(Dvec2 ~ pmat[, 2], family = binomial(link = "probit"), offset = scale(pmat[, 1]),
                control = list(maxit = max_iter))    
    res_pmod <- coef(summary(pmod))
    rownames(res_pmod) <- c("intercept", "trt")
  } else if(opt_std==F & opt_offset){
    pmod <- glm(Dvec2 ~ pmat[, 2], family = binomial(link = "probit"), offset = pmat[, 1],
                control = list(maxit = max_iter))    
    res_pmod <- coef(summary(pmod))
    rownames(res_pmod) <- c("intercept", "trt")
  } else if(opt_std & opt_offset==F){
    pmod <- glm(Dvec2 ~ scale(pmat[,1]) + pmat[, 2], family = binomial(link = "probit"),
                control = list(maxit = max_iter))    
    res_pmod <- coef(summary(pmod))
    rownames(res_pmod) <- c("intercept", "baseline", "trt")
  } else if(opt_std==F & opt_offset==F){
    pmod <- glm(Dvec2 ~ pmat[,1] + pmat[, 2], family = binomial(link = "probit"),
                control = list(maxit = max_iter))    
    res_pmod <- coef(summary(pmod))
    rownames(res_pmod) <- c("intercept", "baseline", "trt")
  }
  #list(pi_hat = pi_hat, res_lmod = res_lmod, sig_hat = sig_hat, res_pmod = res_pmod)
  return(list(pi_hat = pi_hat, res_lmod = res_lmod,
              sig_hat = sig_hat, res_pmod = res_pmod, resid_vec = resid_vec))
}

#############################################################################
###   2. Function for conducting a bootstrap for the first scheme         ###
#############################################################################

boots_fn1 <- function(B, gam0_hat, gam2_hat, beta_intercept, beta_base, beta_hat, delta_hat,
                      scaled_baseline, pmat, Rvec, resid_vec, max_iter = 500){
  
  n <- length(Rvec)
  centered_resid <- resid_vec - mean(resid_vec)
  ind_not_miss <- which(Rvec>=0)
  ind_rd <- which(Rvec==0)
  
  boots_gam0 <- boots_gam2 <- boots_beta <- boots_intercept <- boots_delta <- boots_tp <- rep(NA, B)
  
  for(b in 1:B){
    
    vv <- gam0_hat + scaled_baseline + gam2_hat*pmat[,2] + rnorm(n)
    boots_Dvec <- ifelse(vv < 0, 1, 0)
    boots_Dvec2 <- ifelse(boots_Dvec==1, 0, 1)
    pmod <- glm(boots_Dvec2 ~ pmat[,2], family = binomial(link = "probit"),
                offset = scale(pmat[, 1]), control = list(maxit = max_iter))    
    boots_gam0[b] <- coef(summary(pmod))[1,1]
    boots_gam2[b] <- coef(summary(pmod))[2,1]
    
    boots_y <- beta_intercept + (1+beta_base)*pmat[ind_not_miss,1] + beta_hat*pmat[ind_not_miss,2] +
      delta_hat*ifelse(ind_not_miss %in% ind_rd, 1, 0) + sample(centered_resid, size = length(ind_not_miss), replace = T)  
    
    # the difference between responses would be used as new response
    new_resp <- boots_y - pmat[ind_not_miss, 1]  
    
    # set a new design matrix
    new_dmat <- cbind(ifelse(Rvec[ind_not_miss]==0, 1, 0), pmat[ind_not_miss,])
    lmod <- lm(new_resp ~ new_dmat)
    boots_intercept[b] <- coef(summary(lmod))[1,1]
    boots_delta[b] <- coef(summary(lmod))[2,1]
    boots_beta[b] <- coef(summary(lmod))[4,1]
    
    boots_upp <- boots_gam0[b] + scaled_baseline + boots_gam2[b]
    boots_lwr <- boots_gam0[b] + scaled_baseline
    boots_tp[b] <- boots_beta[b] + boots_delta[b]*mean(pnorm(boots_upp) - pnorm(boots_lwr))
    
  }
  
  return(boots_tp)
}

