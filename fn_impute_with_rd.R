#############################################################################
###                                                                       ###
###   Objective: Imputation methods with retrieved dropout (RD) data      ###
###                                                                       ###
###   Author: Sangyoon Yi (sayi@okstate.edu)                              ###
###           Myeongjong Kang (mkangstat@gmail.com)                       ###
###                                                                       ###
###   Key functions: ancova_rtb, ancova_wo, ancova_rd_simple              ###
###                                                                       ###
###   Reference: Estimation of treatment effect in clinical trials of     ###
###                         continuous endpoints with retrieved dropouts  ###
###                                                                       ###
#############################################################################

library(rstatix)
library(emmeans)

#############################################################################
###   1. Return-to-baseline imputation                                    ###
#############################################################################

# Inputs
# ss: the resulting list object from our sim_model_simple function
# m: the total number of imputed samples to be considered

# Outputs
# estimate: the resulting estimates of \beta_{1}
# std: the standard error of the estimate
# z-value: the test statistic for H_{0}: \beta_{1} = 0
# Pr(>|z|): the corresponding p-value

ancova_rtb <- function(ss, m, sig_lv = 0.05)
{
  Y_baseline    <- ss$pmat[, 1]
  Y_endpoint    <- ss$obs_Y - Y_baseline
  trt           <- ss$pmat[, 2]
  # disc_smed     <- as.numeric(ss$Rvec >= 1) # = 1 (a subject continues study med) or 0 (a subject discontinues study med)
  disc_miss     <- as.numeric(ss$Rvec >= 0) # = 1 (a subject has observed endpoint) or 0 (a subject has missing endpoint)
  
  data_obs      <- data.frame(baseline = Y_baseline[disc_miss == 1], group = trt[disc_miss == 1], endpoint = Y_endpoint[disc_miss == 1])
  data_new      <- data.frame(baseline = Y_baseline[disc_miss == 0], group = trt[disc_miss == 0], endpoint = NA) # NA instead of Y_endpoint[disc_miss == 0]
  fit.out       <- aov(endpoint ~ baseline + group, data = data_obs)
  std.ept       <- sigma(fit.out)
  
  # Repeat the below m times
  estimates     <- rep(NA, m)
  stderrors     <- rep(NA, m)
  for(i in 1:m) {
    
    data_new$endpoint <- rnorm(n = nrow(data_new), mean = rep(0, nrow(data_new)), sd = rep(std.ept, nrow(data_new)))
    ancv.fit      <- summary(lm(endpoint ~ baseline + group, data = rbind(data_obs, data_new)))
    
    estimates[i]  <- ancv.fit$coefficients["group", "Estimate"]
    stderrors[i]  <- ancv.fit$coefficients["group", "Std. Error"]
  }
  
  avg.ests <- mean(estimates)
  std.ests <- sqrt( mean(stderrors^2) + (1 + 1/m) * var(estimates) )
  ts <- avg.ests / std.ests
  lwr_bd <- avg.ests - qnorm(1-sig_lv/2)*std.ests
  upp_bd <- avg.ests + qnorm(1-sig_lv/2)*std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, pval = 2 * (1 - pnorm(q = abs(ts))),
              conf.int = c(lwr_bd, upp_bd)))
}

# ancova_rtb(ss = ss, m = 1000)

#############################################################################
###   2. Washout imputation                                               ###
#############################################################################

# Inputs
# ss: the resulting list object from our sim_model_simple function
# m: the total number of imputed samples to be considered

# Outputs
# estimate: the resulting estimates of \beta_{1}
# std: the standard error of the estimate
# z-value: the test statistic for H_{0}: \beta_{1} = 0
# Pr(>|z|): the corresponding p-value

ancova_wo <- function(ss, m, sig_lv = 0.05)
{
  Y_baseline    <- ss$pmat[, 1]
  Y_endpoint    <- ss$obs_Y - Y_baseline
  trt           <- ss$pmat[, 2]
  disc_smed     <- as.numeric(ss$Rvec >= 1) # = 1 (a subject continues study med) or 0 (a subject discontinues study med)
  # disc_miss     <- as.numeric(ss$Rvec >= 0) # = 1 (a subject has observed endpoint) or 0 (a subject has missing endpoint)
  
  data_obs      <- data.frame(baseline = Y_baseline[disc_smed == 1], group = trt[disc_smed == 1], endpoint = Y_endpoint[disc_smed == 1])
  data_new      <- data.frame(baseline = Y_baseline[disc_smed == 0], group = trt[disc_smed == 0], endpoint = NA) # NA instead of Y_endpoint[disc_smed == 0]
  fit.out       <- aov(endpoint ~ baseline, data = data_obs[data_obs$group == 0, , drop = FALSE]) # Use the data from placebo
  std.ept       <- sigma(fit.out)
  
  # Repeat the below m times
  estimates     <- rep(NA, m)
  stderrors     <- rep(NA, m)
  for(i in 1:m) {
    
    data_new$endpoint <- rnorm(n = nrow(data_new), mean = predict(fit.out, data_new), sd = rep(std.ept, nrow(data_new)))
    ancv.fit      <- summary(lm(endpoint ~ baseline + group, data = rbind(data_obs, data_new)))
    
    estimates[i]  <- ancv.fit$coefficients["group", "Estimate"]
    stderrors[i]  <- ancv.fit$coefficients["group", "Std. Error"]
  }
  
  avg.ests <- mean(estimates)
  std.ests <- sqrt( mean(stderrors^2) + (1 + 1/m) * var(estimates) )
  ts <- avg.ests / std.ests
  lwr_bd <- avg.ests - qnorm(1-sig_lv/2)*std.ests
  upp_bd <- avg.ests + qnorm(1-sig_lv/2)*std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, pval = 2 * (1 - pnorm(q = abs(ts))),
              conf.int = c(lwr_bd, upp_bd)))
}

#############################################################################
###   3. Retrieved-dropout imputation                                     ###
#############################################################################

# Inputs
# ss: the resulting list object from our sim_model_simple function
# m: the total number of imputed samples to be considered

# Outputs
# estimate: the resulting estimates of \beta_{1}
# std: the standard error of the estimate
# z-value: the test statistic for H_{0}: \beta_{1} = 0
# Pr(>|z|): the corresponding p-value

ancova_rd_simple <- function(ss, m, sig_lv = 0.05)
{
  Y_baseline    <- ss$pmat[, 1]
  Y_endpoint    <- ss$obs_Y
  trt           <- ss$pmat[, 2]
  # disc_smed     <- as.numeric(ss$Rvec >= 1) # = 1 (a subject continues study med) or 0 (a subject discontinues study med)
  # disc_miss     <- as.numeric(ss$Rvec >= 0) # = 1 (a subject has observed endpoint) or 0 (a subject has missing endpoint)
  
  data_obs      <- data.frame(baseline = Y_baseline[ss$Rvec == 1], group = trt[ss$Rvec == 1], endpoint = Y_endpoint[ss$Rvec == 1])
  data_rd       <- data.frame(baseline = Y_baseline[ss$Rvec == 0], group = trt[ss$Rvec == 0], endpoint = Y_endpoint[ss$Rvec == 0])
  data_new      <- data.frame(baseline = Y_baseline[ss$Rvec == -1], group = trt[ss$Rvec == -1], endpoint = NA) # NA instead of Y_endpoint[ss$Rvec == -1]
  fit.out       <- aov(endpoint ~ baseline + group, data = data_rd)
  std.ept       <- sigma(fit.out)
  
  # Repeat the below m times
  estimates     <- rep(NA, m)
  stderrors     <- rep(NA, m)
  for(i in 1:m) {
    
    data_new$endpoint <- rnorm(n = nrow(data_new), mean = predict(fit.out, data_new), sd = rep(std.ept, nrow(data_new)))
    ancv.fit      <- summary(lm(endpoint ~ baseline + group, data = rbind(data_obs, data_rd, data_new)))
    
    estimates[i]  <- ancv.fit$coefficients["group", "Estimate"]
    stderrors[i]  <- ancv.fit$coefficients["group", "Std. Error"]
  }
  
  avg.ests <- mean(estimates)
  std.ests <- sqrt( mean(stderrors^2) + (1 + 1/m) * var(estimates) )
  ts <- avg.ests / std.ests
  lwr_bd <- avg.ests - qnorm(1-sig_lv/2)*std.ests
  upp_bd <- avg.ests + qnorm(1-sig_lv/2)*std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, pval = 2 * (1 - pnorm(q = abs(ts))),
              conf.int = c(lwr_bd, upp_bd)))
}
