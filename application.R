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

rm(list = ls())

library(haven)
library(tidyverse)
library(ggplot2)

source("fn_main.R") # load the source code we have

##############################################################################
###
##############################################################################

dat <- read_sas("chapter15_example.sas7bdat")

# let each row be the obsn' for each patient at the endpoint 
df0 <- dat %>%
  group_by(PATIENT) %>%
  slice_max(order_by = VISIT, n = 1, with_ties = FALSE) %>%
  ungroup()

##############################################################################
###
##############################################################################

nrep <- 500 # the total number of runs
k <- 10  # percentage of subjects to select from top (k + 10)%
Rf <- 0.5 # control parameter for persisting drug effect

##############################################################################
###
##############################################################################

res_arr <- array(NA, dim = c(nrep, 2, 4))
for(j in 1:nrep){
  
  set.seed(2025 + j)
  print(paste0("Iteration ", j, " started: ", Sys.time()))

  df <- as.data.frame(df0)
  
  # Filter and rank
  df_top <- df[df$VISIT == 7,] %>%
    filter(THERAPY == "DRUG") %>%
    arrange(change) %>%  # ascending order: best change is smallest
    mutate(rank = row_number()) %>%
    { 
      n_total <- nrow(.)                             # total number of DRUG subjects
      n_top <- ceiling((k + 10) / 100 * n_total)  # top (k+10)% count
      n_k <- ceiling(k / 100 * n_total)              # sample size to draw
      top_k10 <- slice_head(., n = n_top)        # top (k+10)% performers
      slice_sample(top_k10, n = n_k)                 # randomly pick k% from them
    }
  
  df_bottom <- df[df$VISIT == 7,] %>%
    filter(THERAPY == "PLACEBO") %>%
    arrange(desc(change)) %>%  # ascending order: best change is smallest
    mutate(rank = row_number()) %>%
    { 
      n_total <- nrow(.)                             # total number of DRUG subjects
      n_bottom <- ceiling((k + 10) / 100 * n_total)  # top (k+10)% count
      n_k <- ceiling(k / 100 * n_total)              # sample size to draw
      top_k10 <- slice_head(., n = n_bottom)        # top (k+10)% performers
      slice_sample(top_k10, n = n_k)                 # randomly pick k% from them
    }
  
  Rvec <- ifelse(df$VISIT == 7, 1, -1)
  
  Rvec[which(df$PATIENT %in% df_top$PATIENT)] <- 0
  df$change[which(df$PATIENT %in% df_top$PATIENT)] <- df$change[which(df$PATIENT %in% df_top$PATIENT)] * Rf
  df$HAMDTL17[which(df$PATIENT %in% df_top$PATIENT)] <- df$basval[which(df$PATIENT %in% df_top$PATIENT)] + df$change[which(df$PATIENT %in% df_top$PATIENT)]

  Rvec[which(df$PATIENT %in% df_bottom$PATIENT)] <- 0
  df$change[which(df$PATIENT %in% df_bottom$PATIENT)] <- df$change[which(df$PATIENT %in% df_bottom$PATIENT)] * Rf
  df$HAMDTL17[which(df$PATIENT %in% df_top$PATIENT)] <- df$basval[which(df$PATIENT %in% df_top$PATIENT)] + df$change[which(df$PATIENT %in% df_top$PATIENT)]
  
  pmat <- cbind(df$basval, ifelse(df$THERAPY == "DRUG", 1, 0))
  obs_y <- as.vector(df$HAMDTL17)
  obs_y[Rvec < 0] <- 0 
  
  fit <- find_mle(pmat, obs_y, Rvec, max_iter = 500, opt_std = T, opt_offset = T)
  intercept_hat <- fit$res_lmod[1,1] ; delta_hat <- fit$res_lmod[2,1] 
  beta_base <- fit$res_lmod[3,1] ; beta_hat <- fit$res_lmod[4,1]
  gam0_hat <- fit$res_pmod[1,1] ; gam2_hat <- fit$res_pmod[2,1]
  
  Ytilde <- as.vector(scale(pmat[,1]))
  res_arr[j,1,1] <-  beta_hat +  delta_hat*mean(pnorm(gam0_hat + Ytilde + gam2_hat) - pnorm(gam0_hat + Ytilde))
  resid_vec <- fit$resid_vec
  
  # do bootstrap with the 1st scheme
  bb <- boots_fn1(B = 1000, gam0_hat, gam2_hat, intercept_hat, beta_base, beta_hat, delta_hat, 
                  Ytilde, pmat, Rvec, resid_vec, max_iter = 500)
  
  upp_bd <- 2*res_arr[j,1,1] - unname(quantile(bb, prob = 0.05/2))
  lwr_bd <- 2*res_arr[j,1,1] - unname(quantile(bb, prob = 1 - 0.05/2))
  
  if((lwr_bd <= 0) & (upp_bd >= 0)){
    res_arr[j,2,1] <- F
  } else{
    res_arr[j,2,1] <- T
  }
  
  # competing methods
  dat <- list(pmat = pmat, Rvec = Rvec, obs_Y = obs_y)
  
  rtb_obj <- ancova_rtb(dat, m = 1000, sig_lv = 0.05) # return-to-baseline
  res_arr[j,1,2] <- rtb_obj$bhat
  if(rtb_obj$pval < 0.05){
    res_arr[j,2,2] <- T
  } else{
    res_arr[j,2,2] <- F
  }
  
  ws_obj <- ancova_wo(dat, m = 1000, sig_lv = 0.05) # washout
  res_arr[j,1,3] <- ws_obj$bhat
  if(ws_obj$pval < 0.05){
    res_arr[j,2,3] <- T
  } else{
    res_arr[j,2,3] <- F
  }
  
  rd_obj <- ancova_rd_simple(dat, m = 1000, sig_lv = 0.05)  # rd imputation
  res_arr[j,1,4] <- rd_obj$bhat
  if(rd_obj$pval < 0.05){
    res_arr[j,2,4] <- T
  } else{
    res_arr[j,2,4] <- F
  }
  
}

##############################################################################
###
##############################################################################

# save.image(file = "appout.RData")

res_df <- data.frame("Estimate" = c(res_arr[,1,1], 
                                    res_arr[,1,2], 
                                    res_arr[,1,3], 
                                    res_arr[,1,4]),
                     "Method" = c(rep("Our proposed method",nrep), 
                                  rep("RTB imputation",nrep), 
                                  rep("Washout imputation",nrep), 
                                  rep("RD imputation", nrep)))

p01 <- res_df %>% 
  mutate(Method = factor(Method, levels = c("Our proposed method", 
                                            "RD imputation", 
                                            "RTB imputation", 
                                            "Washout imputation"))) %>% 
  ggplot(aes(x = Method, y = Estimate, fill = Method)) + 
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Method") + ylab("Treatment effect estimate under the TP strategy") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

ggsave(filename = "appout_boxplot_Rf=0.5.pdf", plot = p01, width = 10, height = 4) # Change Rf (from 0.5 to 0) to obtain appout_boxplot_Rf=0.pdf