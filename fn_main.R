#############################################################################
###                                                                       ###
###   Objective: Running simulations and summarizing results              ###
###                                                                       ###
###   Author: Sangyoon Yi (sayi@okstate.edu)                              ###
###           Myeongjong Kang (mkangstat@gmail.com)                       ###
###                                                                       ###
###   Key functions: sim_fn                                               ###
###                                                                       ###
###   Reference: Estimation of treatment effect in clinical trials of     ###
###                         continuous endpoints with retrieved dropouts  ###
###                                                                       ###
#############################################################################

source("fn_gen_simdata.R")
source("fn_ancova_with_rd.R")
source("fn_impute_with_rd.R")

#############################################################################
###   1. function to simulate a single time                               ###
#############################################################################

sim_fn <- function(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob,
                   gam0, gam1 = 1, gam2, seednum = NULL, sig_lv = 0.05){
  
  #beta2 <- 0
  if(!is.null(seednum)){
    set.seed(seednum)
  }
  
  bvec <- c(1, beta2)
  gvec <- c(gam0, gam1, gam2)
  delta <- kappa * beta2
  true_tp <- bvec[2] + delta*(pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
  
  neg_indc <- T
  
  while(neg_indc){
    dat <- sim_model(n, ntrt, qprob, delta, bvec, gvec, mu_base, sd_base, sd_eps, sd_u = 1)  
    if(min(dat$obs_Y) >= 0){
      neg_indc <- F
    }
  }
  
  res_mat <- matrix(NA, nrow = 4, ncol = 4)  
  rownames(res_mat) <- c("new", "rtb", "wash", "rd_imp")
  colnames(res_mat) <- c("bhat", "covered", "length", "rejected")
  
  # 1. Our method
  fit <- find_mle(pmat = dat$pmat, obs_y = dat$obs_Y, Rvec = dat$Rvec, max_iter = 500, opt_std = T, opt_offset = T)
  intercept_hat <- fit$res_lmod[1,1] ; delta_hat <- fit$res_lmod[2,1]
  beta_base <- fit$res_lmod[3,1] ; beta_hat <- fit$res_lmod[4,1]
  gam0_hat <- fit$res_pmod[1,1] ; gam2_hat <- fit$res_pmod[2,1]
  
  Ytilde <- as.vector(scale(dat$pmat[,1]))
  #upp <- gam0_hat + Ytilde + gam2_hat
  #lwr <- gam0_hat + Ytilde
  res_mat[1,1] <-  beta_hat +  delta_hat*mean(pnorm(gam0_hat + Ytilde + gam2_hat) - pnorm(gam0_hat + Ytilde))
  #res_mat[1,1] <- bhat_tp
  
  # do bootstrap with the 1st scheme
  bb <- boots_fn1(B = 1000, gam0_hat, gam2_hat, intercept_hat, beta_base, beta_hat, delta_hat,
                  Ytilde, dat$pmat, dat$Rvec, fit$resid_vec, max_iter = 500)
  
  upp_bd <- 2*res_mat[1,1] - unname(quantile(bb, prob = sig_lv/2))
  lwr_bd <- 2*res_mat[1,1] - unname(quantile(bb, prob = 1 - sig_lv/2))
  #c(lwr_bd, upp_bd)
  if(lwr_bd<=true_tp && upp_bd>=true_tp){
    res_mat[1,2] <- 1
  } else{
    res_mat[1,2] <- 0
  }
  res_mat[1,3] <- upp_bd - lwr_bd
  boots.pval <- (sum(abs(bb - mean(bb)) >= abs(res_mat[1,1])) + 1)/(1000 + 1)
  res_mat[1,4] <- ifelse(boots.pval<=sig_lv, 1, 0)
  
  # 2. competing methods
  rtb_obj <- ancova_rtb(dat, m = 1000, sig_lv) # return-to-baseline
  res_mat[2,1] <- rtb_obj$bhat
  if(rtb_obj$conf.int[1]<=true_tp && rtb_obj$conf.int[2]>=true_tp){
    res_mat[2,2] <- 1
  } else{
    res_mat[2,2] <- 0
  }
  res_mat[2,3] <- diff(rtb_obj$conf.int)
  res_mat[2,4] <- ifelse(rtb_obj$pval<=sig_lv, 1, 0)
  
  ws_obj <- ancova_wo(dat, m = 1000, sig_lv) # washout
  res_mat[3,1] <- ws_obj$bhat
  if(ws_obj$conf.int[1]<=true_tp && ws_obj$conf.int[2]>=true_tp){
    res_mat[3,2] <- 1
  } else{
    res_mat[3,2] <- 0
  }
  res_mat[3,3] <- diff(ws_obj$conf.int)
  res_mat[3,4] <- ifelse(ws_obj$pval<=sig_lv, 1, 0)
  
  rd_obj <- ancova_rd_simple(dat, m = 1000, sig_lv)  # rd imputation
  res_mat[4,1] <- rd_obj$bhat
  if(rd_obj$conf.int[1]<=true_tp && rd_obj$conf.int[2]>=true_tp){
    res_mat[4,2] <- 1
  } else{
    res_mat[4,2] <- 0
  }  
  res_mat[4,3] <- diff(rd_obj$conf.int)
  res_mat[4,4] <- ifelse(rd_obj$pval<=sig_lv, 1, 0)
  # res_mat ; true_tp
  #return(list(pi_hat = fit$pi_hat, sig_hat = fit$sig_hat,
  #            delta_hat = fit$res_lmod[1,1], pval_delta = fit$res_lmod[1,4],
  #            pmod = fit$res_pmod, res_mat = res_mat))
  return(res_mat)
}