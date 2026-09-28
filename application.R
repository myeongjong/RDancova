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

source("fn_main.R")

##############################################################################
### Read the original longitudinal dataset                                 ###
##############################################################################

dat_long  <- read_sas("chapter15_example.sas7bdat")

### Retain the last available record for each subject
df0       <- dat_long %>% group_by(PATIENT) %>% slice_max(order_by = VISIT, n = 1, with_ties = FALSE) %>% ungroup()

##############################################################################
### Simulation parameters                                                  ###
##############################################################################

nrep      <- 500
k         <- 10
Rf        <- 0.50
Bf        <- 0.00

### Rf: proportion of the original treatment effect retained
### Bf: proportion of the rescue-treatment benefit received
if(length(Rf) != 1 || !is.finite(Rf) || Rf < 0 || Rf > 1) stop("Rf must be a single value between 0 and 1.")
if(length(Bf) != 1 || !is.finite(Bf) || Bf < 0 || Bf > 1) stop("Bf must be a single value between 0 and 1.")

##############################################################################
### Estimate the treatment effect and residual SD from completers (Wang et al. use an MMRM; this ANCOVA is a final-visit approximation)
##############################################################################

df_comp   <- as.data.frame(df0)
df_comp   <- df_comp[df_comp$VISIT == 7, , drop = FALSE]
df_comp$trt <- ifelse(df_comp$THERAPY == "DRUG", 1, 0)

fit_rep   <- lm(change ~ basval + trt, data = df_comp)
EC        <- unname(coef(fit_rep)["trt"])
sigma_rep <- summary(fit_rep)$sigma

if(!is.finite(EC)) stop("The completer-based treatment effect could not be estimated.")
if(!is.finite(sigma_rep) || sigma_rep <= 0) stop("The completer-based residual standard deviation must be positive.")

##############################################################################
### Main                                                                   ###
##############################################################################

### Store treatment-effect estimates and rejection indicators
res_arr   <- array(NA, dim = c(nrep, 2, 4))

for (j in 1:nrep) {
  
  set.seed(2025 + j)
  print(paste0("Iteration ", j, " started: ", Sys.time()))
  
  df        <- as.data.frame(df0)
  
  ### Select top-performing DRUG subjects
  df_top    <- df[df$VISIT == 7, ] %>% filter(THERAPY == "DRUG") %>% arrange(change) %>% mutate(rank = row_number()) %>%
    {
      n_total     <- nrow(.)
      n_top       <- ceiling((k + 10) / 100 * n_total)
      n_k         <- ceiling(k / 100 * n_total)
      top_k10     <- slice_head(., n = n_top)
      slice_sample(top_k10, n = n_k)
    }
  
  ### Select bottom-performing PLACEBO subjects
  df_bottom <- df[df$VISIT == 7, ] %>% filter(THERAPY == "PLACEBO") %>% arrange(desc(change)) %>% mutate(rank = row_number()) %>%
    {
      n_total     <- nrow(.)
      n_bottom    <- ceiling((k + 10) / 100 * n_total)
      n_k         <- ceiling(k / 100 * n_total)
      bottom_k10  <- slice_head(., n = n_bottom)
      slice_sample(bottom_k10, n = n_k)
    }
  
  ### Rvec = 1: completer
  ### Rvec = 0: retrieved dropout
  ### Rvec = -1: non-retrieved dropout with a missing endpoint
  Rvec      <- ifelse(df$VISIT == 7, 1, -1)
  
  ### Wang Setting 1 for selected experimental subjects
  idx_top   <- which(df$PATIENT %in% df_top$PATIENT)
  mu_top    <- df$HAMDTL17[idx_top] + abs(df$HAMDTL17[idx_top] - df$basval[idx_top]) * (1 - Rf)
  
  df$HAMDTL17[idx_top]  <- rnorm(length(idx_top), mean = mu_top, sd = sigma_rep)
  df$change[idx_top]    <- df$HAMDTL17[idx_top] - df$basval[idx_top]
  
  Rvec[idx_top] <- 0
  
  ### Wang Setting 2 for selected placebo subjects
  idx_bottom  <- which(df$PATIENT %in% df_bottom$PATIENT)
  mu_bottom   <- df$HAMDTL17[idx_bottom] - abs(EC) * Bf
  
  df$HAMDTL17[idx_bottom] <- rnorm(length(idx_bottom), mean = mu_bottom, sd = sigma_rep)
  df$change[idx_bottom]   <- df$HAMDTL17[idx_bottom] - df$basval[idx_bottom]
  
  Rvec[idx_bottom] <- 0
  
  ### Construct the analysis inputs
  pmat  <- cbind(df$basval, ifelse(df$THERAPY == "DRUG", 1, 0))
  obs_y <- as.vector(df$HAMDTL17)
  
  ### Zero is used only as a placeholder for a missing endpoint
  obs_y[Rvec == -1] <- 0
  
  ss_obj    <- list(pmat = pmat, Rvec = Rvec, obs_Y = obs_y)
  
  ### Proposed method
  fit       <- find_mle(pmat, obs_y, Rvec, max_iter = 500, opt_std = TRUE, opt_offset = TRUE)
  
  intercept_hat <- fit$res_lmod[1, 1]
  delta_hat     <- fit$res_lmod[2, 1]
  beta_base     <- fit$res_lmod[3, 1]
  beta_hat      <- fit$res_lmod[4, 1]
  gam0_hat      <- fit$res_pmod[1, 1]
  gam2_hat      <- fit$res_pmod[2, 1]
  
  Ytilde    <- as.vector(scale(pmat[, 1]))
  
  res_arr[j, 1, 1] <- beta_hat + delta_hat * mean(pnorm(gam0_hat + Ytilde + gam2_hat) - pnorm(gam0_hat + Ytilde))
  
  resid_vec <- fit$resid_vec
  
  ### Bootstrap inference for the proposed method
  bb      <- boots_fn1(B = 1000, gam0_hat, gam2_hat, intercept_hat, beta_base, beta_hat, delta_hat, Ytilde, pmat, Rvec, resid_vec, max_iter = 500)
  upp_bd  <- 2 * res_arr[j, 1, 1] - unname(quantile(bb, prob = 0.05/2))
  lwr_bd  <- 2 * res_arr[j, 1, 1] - unname(quantile(bb, prob = 1 - 0.05/2))
  
  res_arr[j, 2, 1] <- !(lwr_bd <= 0 && upp_bd >= 0)
  
  ### Return-to-baseline imputation
  rtb_obj <- ancova_rtb(ss = ss_obj, m = 1000, sig_lv = 0.05)
  
  res_arr[j, 1, 2] <- rtb_obj$bhat
  res_arr[j, 2, 2] <- rtb_obj$pval < 0.05
  
  ### Washout imputation
  ws_obj  <- ancova_wo(ss = ss_obj, m = 1000, sig_lv = 0.05)
  
  res_arr[j, 1, 3] <- ws_obj$bhat
  res_arr[j, 2, 3] <- ws_obj$pval < 0.05
  
  ### Retrieved-dropout imputation
  rd_obj  <- ancova_rd_simple(ss = ss_obj, m = 1000, sig_lv = 0.05)
  
  res_arr[j, 1, 4] <- rd_obj$bhat
  res_arr[j, 2, 4] <- rd_obj$pval < 0.05
}

##############################################################################
### Create the boxplot                                                     ###
##############################################################################

res_df <- data.frame(Estimate = c(res_arr[, 1, 1], 
                                  res_arr[, 1, 2], 
                                  res_arr[, 1, 3], 
                                  res_arr[, 1, 4]), 
                     Method = c(rep("Our proposed method", nrep), 
                                rep("RTB imputation", nrep), 
                                rep("Washout imputation", nrep), 
                                rep("RD imputation", nrep)))

res_df$Method <- factor(res_df$Method, levels = c("Our proposed method", "RD imputation", "RTB imputation", "Washout imputation"))

p01 <- ggplot(res_df, aes(x = Method, y = Estimate, fill = Method)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Method") +
  ylab("Treatment effect estimate under the TP strategy") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12))

##############################################################################
### Save the result and boxplot                                            ###
##############################################################################

result_file <- sprintf("appout_Rf=%.2f_Bf=%.2f.RData", Rf, Bf)
figure_file <- sprintf("appout_main_boxplot_Rf=%.2f_Bf=%.2f.pdf", Rf, Bf)

save(res_df, res_arr, file = result_file)
ggsave(filename = figure_file, plot = p01, width = 10, height = 4)

##############################################################################
### Create boxplots for the supplementary material                         ###
##############################################################################

# rm(list = ls())
# 
# library(haven)
# library(tidyverse)
# library(ggplot2)
# 
# Rf        <- 1.00
# Bf        <- 0.00
# 
# result_file <- sprintf("appout_Rf=%.2f_Bf=%.2f.RData", Rf, Bf)
# figure_file <- sprintf("appout_main_boxplot_Rf=%.2f_Bf=%.2f_V2.pdf", Rf, Bf)
# 
# load(file = result_file)
# 
# if(Rf < 1) {
#   
#   p01 <- ggplot(res_df, aes(x = Method, y = Estimate, fill = Method)) +
#     geom_boxplot() +
#     scale_y_continuous(breaks = seq(-12, 4, by = 1)) +
#     coord_flip(ylim = c(-12, 4)) +
#     scale_fill_brewer(palette = "Dark2") +
#     xlab("Method") +
#     ylab("") +
#     theme_bw() +
#     theme(legend.position = "none", legend.title = element_blank(),
#           axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12),
#           axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12))
#   
# } else {
#   
#   p01 <- ggplot(res_df, aes(x = Method, y = Estimate, fill = Method)) +
#     geom_boxplot() +
#     scale_y_continuous(breaks = seq(-12, 4, by = 1)) +
#     coord_flip(ylim = c(-12, 4)) +
#     scale_fill_brewer(palette = "Dark2") +
#     xlab("Method") +
#     ylab("Treatment effect estimate under the TP strategy") +
#     theme_bw() +
#     theme(legend.position = "none", legend.title = element_blank(),
#           axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12),
#           axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12))
# }
# 
# ggsave(filename = figure_file, plot = p01, width = 10, height = 2)