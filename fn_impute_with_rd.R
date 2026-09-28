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

#############################################################################
###   1. Return-to-baseline imputation                                    ###
#############################################################################

### Inputs
###   ss: List returned by sim_model_simple(), containing baseline values, treatment groups, observed endpoints, and endpoint-observation status.
###   m: Number of imputed datasets to generate.
###   sig_lv: Significance level for the two-sided test and confidence interval (default: 0.05, giving a 95% confidence interval).
###
### Output
###   A list containing:
###     bhat: Pooled ANCOVA treatment-effect estimate.
###     std.err: Multiple-imputation standard error.
###     ts: Test statistic for a treatment effect of zero.
###     df: Degrees of freedom used for inference.
###     pval: Two-sided p-value.
###     conf.int: Lower and upper confidence limits.

ancova_rtb <- function(ss, m, sig_lv = 0.05)
{
  data_all      <- data.frame(baseline = ss$pmat[, 1], group = ss$pmat[, 2], endpoint = ss$obs_Y, Rvec = ss$Rvec)
  miss          <- data_all$Rvec == -1
  
  mean_baseline <- mean(data_all$baseline)
  
  estimates     <- rep(NA, m)
  stderrors     <- rep(NA, m)
  for (i in 1:m) {
    
    data_imp      <- data_all
    endpoint_mar  <- data_imp$endpoint
    
    for (g in c(0, 1)) {
      
      data_obs_g    <- data_all[data_all$group == g & data_all$Rvec >= 0, , drop = FALSE]
      mis_g         <- which(data_all$group == g & data_all$Rvec == -1)
      
      if (length(mis_g) > 0) {
        
        fit.out       <- lm(endpoint ~ baseline, data = data_obs_g)
        X             <- model.matrix(fit.out)
        beta_hat      <- coef(fit.out)
        SSE           <- sum(residuals(fit.out)^2)
        df_res        <- df.residual(fit.out)
        
        if (fit.out$rank < ncol(X) || df_res <= 0 || !is.finite(SSE) || SSE <= 0) stop(paste0("RTB MAR model cannot be estimated in treatment group ", g, "."))
        
        sigma2_draw <- SSE / rchisq(1, df = df_res)
        XtX_inv     <- solve(crossprod(X))
        cov_beta    <- sigma2_draw * XtX_inv
        cov_beta    <- (cov_beta + t(cov_beta)) / 2
        
        if (!all(is.finite(cov_beta)) || min(diag(cov_beta)) <= 0) stop(paste0("Invalid RTB MAR coefficient covariance in treatment group ", g, "."))
        
        chol_beta   <- tryCatch(chol(cov_beta), error = function(e) NULL)
        if (is.null(chol_beta)) {
          jitter      <- 1e-10 * max(diag(cov_beta))
          chol_beta   <- tryCatch(chol(cov_beta + diag(jitter, ncol(cov_beta))), error = function(e) stop(paste0("RTB MAR covariance remains non-positive-definite in treatment group ", g, ".")))
        }
        
        beta_draw   <- beta_hat + as.vector(t(chol_beta) %*% rnorm(length(beta_hat)))
        
        X_new       <- model.matrix(~ baseline, data = data_all[mis_g, , drop = FALSE])
        mu_draw     <- as.vector(X_new %*% beta_draw)
        
        endpoint_mar[mis_g] <- rnorm(length(mis_g), mean = mu_draw, sd = sqrt(sigma2_draw))
      }
    }
    
    endpoint_rtb <- endpoint_mar
    for (g in c(0, 1)) {
      
      idx_g         <- data_all$group == g
      mis_g         <- which(data_all$group == g & data_all$Rvec == -1)
      
      if (length(mis_g) > 0) {
        
        mean_endpoint_g     <- mean(endpoint_mar[idx_g])
        endpoint_rtb[mis_g] <- endpoint_mar[mis_g] - mean_endpoint_g + mean_baseline
      }
    }
    
    data_imp$endpoint[miss] <- endpoint_rtb[miss]
    
    data_imp$change <- data_imp$endpoint - data_imp$baseline
    ancv.fit        <- lm(change ~ baseline + group, data = data_imp)
    
    estimates[i]    <- coef(ancv.fit)["group"]
    stderrors[i]    <- sqrt(vcov(ancv.fit)["group", "group"])
  }
  
  avg.ests      <- mean(estimates)
  within.var    <- mean(stderrors^2)
  between.var   <- ifelse(m > 1, var(estimates), 0)
  total.var     <- within.var + (1 + 1/m) * between.var
  std.ests      <- sqrt(total.var)
  ts            <- avg.ests / std.ests
  
  if (between.var > 0 && within.var > 0 && m > 1) {
    
    r             <- (1 + 1/m) * between.var / within.var
    df.mi         <- (m - 1) * (1 + 1/r)^2
    crit          <- qt(1 - sig_lv/2, df = df.mi)
    pval          <- 2 * pt(-abs(ts), df = df.mi)
    
  } else {
    
    df.mi         <- Inf
    crit          <- qnorm(1 - sig_lv/2)
    pval          <- 2 * pnorm(-abs(ts))
  }
  
  lwr_bd        <- avg.ests - crit * std.ests
  upp_bd        <- avg.ests + crit * std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, df = df.mi, pval = pval, conf.int = c(lwr_bd, upp_bd)))
}

#############################################################################
###   2. Washout imputation                                               ###
#############################################################################

### Inputs
###   ss: List returned by sim_model_simple(), containing baseline values, treatment groups, observed endpoints, and endpoint-observation status.
###   m: Number of imputed datasets to generate.
###   sig_lv: Significance level for the two-sided test and confidence interval (default: 0.05, giving a 95% confidence interval).
###
### Output
###   A list containing:
###     bhat: Pooled ANCOVA treatment-effect estimate.
###     std.err: Multiple-imputation standard error.
###     ts: Test statistic for a treatment effect of zero.
###     df: Degrees of freedom used for inference.
###     pval: Two-sided p-value.
###     conf.int: Lower and upper confidence limits.

ancova_wo <- function(ss, m, sig_lv = 0.05)
{
  data_all  <- data.frame(baseline = ss$pmat[, 1], group = ss$pmat[, 2], endpoint = ss$obs_Y, Rvec = ss$Rvec)
  miss      <- data_all$Rvec == -1
  
  data_pbo  <- data_all[data_all$group == 0 & data_all$Rvec >= 0, , drop = FALSE]
  
  estimates <- rep(NA, m)
  stderrors <- rep(NA, m)
  for (i in 1:m) {
    
    data_imp <- data_all
    
    if (any(miss)) {
      
      fit.out     <- lm(endpoint ~ baseline, data = data_pbo)
      X           <- model.matrix(fit.out)
      beta_hat    <- coef(fit.out)
      SSE         <- sum(residuals(fit.out)^2)
      df_res      <- df.residual(fit.out)
      
      if (fit.out$rank < ncol(X) || df_res <= 0 || !is.finite(SSE) || SSE <= 0) stop("Washout imputation model cannot be estimated.")
      
      sigma2_draw <- SSE / rchisq(1, df = df_res)
      XtX_inv     <- solve(crossprod(X))
      cov_beta    <- sigma2_draw * XtX_inv
      cov_beta    <- (cov_beta + t(cov_beta)) / 2
      
      if (!all(is.finite(cov_beta)) || min(diag(cov_beta)) <= 0) stop("Invalid washout coefficient covariance.")
      
      chol_beta   <- tryCatch(chol(cov_beta), error = function(e) NULL)
      if (is.null(chol_beta)) {
        
        jitter      <- 1e-10 * max(diag(cov_beta))
        chol_beta   <- tryCatch(chol(cov_beta + diag(jitter, ncol(cov_beta))), error = function(e) stop("Washout covariance remains non-positive-definite."))
      }
      
      beta_draw   <- beta_hat + as.vector(t(chol_beta) %*% rnorm(length(beta_hat)))
      
      X_new       <- model.matrix(~ baseline, data = data_all[miss, , drop = FALSE])
      mu_draw     <- as.vector(X_new %*% beta_draw)
      
      data_imp$endpoint[miss] <- rnorm(sum(miss), mean = mu_draw, sd = sqrt(sigma2_draw))
    }
    
    data_imp$change <- data_imp$endpoint - data_imp$baseline
    ancv.fit        <- lm(change ~ baseline + group, data = data_imp)
    
    estimates[i]    <- coef(ancv.fit)["group"]
    stderrors[i]    <- sqrt(vcov(ancv.fit)["group", "group"])
  }
  
  avg.ests    <- mean(estimates)
  within.var  <- mean(stderrors^2)
  between.var <- ifelse(m > 1, var(estimates), 0)
  total.var   <- within.var + (1 + 1/m) * between.var
  std.ests    <- sqrt(total.var)
  ts          <- avg.ests / std.ests
  
  if (between.var > 0 && within.var > 0 && m > 1) {
    
    r       <- (1 + 1/m) * between.var / within.var
    df.mi   <- (m - 1) * (1 + 1/r)^2
    crit    <- qt(1 - sig_lv/2, df = df.mi)
    pval    <- 2 * pt(-abs(ts), df = df.mi)
    
  } else {
    
    df.mi   <- Inf
    crit    <- qnorm(1 - sig_lv/2)
    pval    <- 2 * pnorm(-abs(ts))
  }
  
  lwr_bd <- avg.ests - crit * std.ests
  upp_bd <- avg.ests + crit * std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, df = df.mi, pval = pval, conf.int = c(lwr_bd, upp_bd)))
}

#############################################################################
###   3. Retrieved-dropout imputation                                     ###
#############################################################################

### Inputs
###   ss: List returned by sim_model_simple(), containing baseline values, treatment groups, observed endpoints, and endpoint-observation status.
###   m: Number of imputed datasets to generate.
###   sig_lv: Significance level for the two-sided test and confidence interval (default: 0.05, giving a 95% confidence interval).
###
### Output
###   A list containing:
###     bhat: Pooled ANCOVA treatment-effect estimate.
###     std.err: Multiple-imputation standard error.
###     ts: Test statistic for a treatment effect of zero.
###     df: Degrees of freedom used for inference.
###     pval: Two-sided p-value.
###     conf.int: Lower and upper confidence limits.

ancova_rd_simple <- function(ss, m, sig_lv = 0.05)
{
  data_all  <- data.frame(baseline = ss$pmat[, 1], group = ss$pmat[, 2], endpoint = ss$obs_Y, Rvec = ss$Rvec)
  
  estimates <- rep(NA, m)
  stderrors <- rep(NA, m)
  for (i in 1:m) {
    
    data_imp <- data_all
    
    for (g in c(0, 1)) {
      
      data_rd_g <- data_all[data_all$Rvec == 0 & data_all$group == g, , drop = FALSE]
      mis_g     <- which(data_all$Rvec == -1 & data_all$group == g)
      if (length(mis_g) > 0) {
        
        fit.out     <- lm(endpoint ~ baseline, data = data_rd_g)
        X           <- model.matrix(fit.out)
        beta_hat    <- coef(fit.out)
        SSE         <- sum(residuals(fit.out)^2)
        df_res      <- df.residual(fit.out)
        
        if (fit.out$rank < ncol(X) || df_res <= 0 || !is.finite(SSE) || SSE <= 0) stop(paste0("RD imputation model cannot be estimated in treatment group ", g, "."))
        
        sigma2_draw <- SSE / rchisq(1, df = df_res)
        XtX_inv     <- solve(crossprod(X))
        cov_beta    <- sigma2_draw * XtX_inv
        cov_beta    <- (cov_beta + t(cov_beta)) / 2
        
        if (!all(is.finite(cov_beta)) || min(diag(cov_beta)) <= 0) stop(paste0("Invalid RD coefficient covariance in treatment group ", g, "."))
        
        chol_beta   <- tryCatch(chol(cov_beta), error = function(e) NULL)
        if (is.null(chol_beta)) {
          
          jitter      <- 1e-10 * max(diag(cov_beta))
          chol_beta   <- tryCatch(chol(cov_beta + diag(jitter, ncol(cov_beta))), error = function(e) stop(paste0("RD covariance remains non-positive-definite in treatment group ", g, ".")))
        }
        
        beta_draw   <- beta_hat + as.vector(t(chol_beta) %*% rnorm(length(beta_hat)))
        
        X_new       <- model.matrix(~ baseline, data = data_all[mis_g, , drop = FALSE])
        mu_draw     <- as.vector(X_new %*% beta_draw)
        
        data_imp$endpoint[mis_g] <- rnorm(length(mis_g), mean = mu_draw, sd = sqrt(sigma2_draw))
      }
    }
    
    data_imp$change <- data_imp$endpoint - data_imp$baseline
    ancv.fit        <- lm(change ~ baseline + group, data = data_imp)
    
    estimates[i]    <- coef(ancv.fit)["group"]
    stderrors[i]    <- sqrt(vcov(ancv.fit)["group", "group"])
  }
  
  avg.ests    <- mean(estimates)
  within.var  <- mean(stderrors^2)
  between.var <- ifelse(m > 1, var(estimates), 0)
  total.var   <- within.var + (1 + 1/m) * between.var
  std.ests    <- sqrt(total.var)
  ts          <- avg.ests / std.ests
  
  if (between.var > 0 && within.var > 0 && m > 1) {
    
    r       <- (1 + 1/m) * between.var / within.var
    df.mi   <- (m - 1) * (1 + 1/r)^2
    crit    <- qt(1 - sig_lv/2, df = df.mi)
    pval    <- 2 * pt(-abs(ts), df = df.mi)
    
  } else {
    
    df.mi   <- Inf
    crit    <- qnorm(1 - sig_lv/2)
    pval    <- 2 * pnorm(-abs(ts))
  }
  
  lwr_bd <- avg.ests - crit * std.ests
  upp_bd <- avg.ests + crit * std.ests
  
  return(list(bhat = avg.ests, std.err = std.ests, ts = ts, df = df.mi, pval = pval, conf.int = c(lwr_bd, upp_bd)))
}