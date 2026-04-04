#############################################################################
###                                                                       ###
###   Objective: Running simulations in parallel on a local machine       ###
###                                                                       ###
###   Author: Sangyoon Yi (sayi@okstate.edu)                              ###
###           Myeongjong Kang (mkangstat@gmail.com)                       ###
###                                                                       ###
###   Key functions: sim_fn, outtab                                       ###
###                                                                       ###
###   Reference: Estimation of treatment effect in clinical trials of     ###
###                         continuous endpoints with retrieved dropouts  ###
###                                                                       ###
#############################################################################

rm(list = ls())

library(doParallel) # load this package for parallel computing
source("fn_main.R") # load the source code we have

##############################################################################
###
##############################################################################

my_seed <- 04012025 # set the number for set.seed()
num_sim <- 5000 # the total number of independent simulation runs

##############################################################################
###
##############################################################################

set.seed(my_seed)
randseeds <- sort(sample(seq(100), size = 12)) * 10000 # sort(sample(x = seq(10000), size = 12))

registerDoParallel(cl <- makeCluster(parallel::detectCores() - 2))

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy01 <- beta2
beta_tp01 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res01 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[1] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy02 <- beta2
beta_tp02 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res02 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[2] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy03 <- beta2
beta_tp03 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res03 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[3] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy04 <- beta2
beta_tp04 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res04 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[4] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy05 <- beta2
beta_tp05 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res05 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[5] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy06 <- beta2
beta_tp06 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res06 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[6] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(-10, -0.5, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy07 <- beta2
beta_tp07 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res07 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[7] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(-10, -1, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy08 <- beta2
beta_tp08 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res08 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[8] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(0, -1, 20, 0.5, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy09 <- beta2
beta_tp09 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res09 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[9] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(-10, -0.5, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy10 <- beta2
beta_tp10 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res10 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[10] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(-10, -1, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy11 <- beta2
beta_tp11 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res11 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[11] + i), error = function(e) NA) }

n <- 100 ; ntrt <- 50 ; args = c(0, -1, 20, 0.5, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy12 <- beta2
beta_tp12 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res12 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[12] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.3, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy13 <- beta2
beta_tp13 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res13 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[1] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.3, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy14 <- beta2
beta_tp14 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res14 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[2] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.3, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy15 <- beta2
beta_tp15 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res15 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[3] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.3, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy16 <- beta2
beta_tp16 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res16 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[4] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.3, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy17 <- beta2
beta_tp17 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res17 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[5] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.3, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy18 <- beta2
beta_tp18 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res18 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[6] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.1, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy19 <- beta2
beta_tp19 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res19 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[1] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.1, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy20 <- beta2
beta_tp20 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res20 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[2] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.1, -0.75, 0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy21 <- beta2
beta_tp21 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res21 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[3] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -0.5, 20, 0.1, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy22 <- beta2
beta_tp22 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res22 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[4] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(-10, -1, 20, 0.1, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy23 <- beta2
beta_tp23 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res23 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[5] + i), error = function(e) NA) }

n <- 200 ; ntrt <- 100 ; args = c(0, -1, 20, 0.1, -0.75, -0.25) ; beta2 = as.numeric(args[1]) ; kappa = as.numeric(args[2]) ; sd_eps = as.numeric(args[3]) ; qprob = as.numeric(args[4]) ; gam0 = as.numeric(args[5]) ; gam2 = as.numeric(args[6]) ; mu_base <- 180 ; sd_base <- 20; gam1 <- 1
beta_hy24 <- beta2
beta_tp24 <- beta2 + (kappa * beta2) * (pnorm((gam0 + gam2)/sqrt(2)) - pnorm((gam0)/sqrt(2)))
res24 <- foreach(i = 1:num_sim) %dopar% { tryCatch(sim_fn(n, ntrt, mu_base, sd_base, beta2, kappa, sd_eps, qprob, gam0, gam1 = gam1, gam2, seednum = randseeds[6] + i), error = function(e) NA) }

stopCluster(cl)

##############################################################################
###
##############################################################################

outtab <- function(res, beta_hy, beta_tp)
{
  print(paste0("Simulation success rate: ", round(sum(!is.na(res)) / length(res) * 100, 1), "% out of ", length(res), " simulations were successful."))
  
  res <- res[!is.na(res)]
  
  bias <- matrix(0, nrow = length(res), ncol = 4) 
  colnames(bias) <- c("new", "rtb", "wash", "rd_imp")
  for(i in 1:nrow(bias)) {
    
    bias[i, "new"] <- bias[i, "new"] + (res[[i]]["new", "bhat"] - beta_tp)
    bias[i, "rtb"] <- bias[i, "rtb"] + (res[[i]]["rtb", "bhat"] - beta_tp)
    bias[i, "wash"] <- bias[i, "wash"] + (res[[i]]["wash", "bhat"] - beta_tp)
    bias[i, "rd_imp"] <- bias[i, "rd_imp"] + (res[[i]]["rd_imp", "bhat"] - beta_tp)
  }
  
  bias <- colMeans(bias) # ; round(bias, 2)
  
  rmse <- matrix(0, nrow = length(res), ncol = 4) 
  colnames(rmse) <- c("new", "rtb", "wash", "rd_imp")
  for(i in 1:nrow(rmse)) {
    
    rmse[i, "new"] <- rmse[i, "new"] + (res[[i]]["new", "bhat"] - beta_tp)^2
    rmse[i, "rtb"] <- rmse[i, "rtb"] + (res[[i]]["rtb", "bhat"] - beta_tp)^2
    rmse[i, "wash"] <- rmse[i, "wash"] + (res[[i]]["wash", "bhat"] - beta_tp)^2
    rmse[i, "rd_imp"] <- rmse[i, "rd_imp"] + (res[[i]]["rd_imp", "bhat"] - beta_tp)^2
  }
  
  rmse <- sqrt(colMeans(rmse)) # ; round(rmse, 2)
  
  covrate <- matrix(0, nrow = length(res), ncol = 4) 
  colnames(covrate) <- c("new", "rtb", "wash", "rd_imp")
  for(i in 1:nrow(covrate)) {
    
    covrate[i, "new"] <- res[[i]]["new", "covered"]
    covrate[i, "rtb"] <- res[[i]]["rtb", "covered"]
    covrate[i, "wash"] <- res[[i]]["wash", "covered"]
    covrate[i, "rd_imp"] <- res[[i]]["rd_imp", "covered"]
  }
  
  covrate <- colMeans(covrate) * 100 # ; round(covrate, 2)
  
  cilength <- matrix(0, nrow = length(res), ncol = 4) 
  colnames(cilength) <- c("new", "rtb", "wash", "rd_imp")
  for(i in 1:nrow(cilength)) {
    
    cilength[i, "new"] <- res[[i]]["new", "length"] 
    cilength[i, "rtb"] <- res[[i]]["rtb", "length"] 
    cilength[i, "wash"] <- res[[i]]["wash", "length"] 
    cilength[i, "rd_imp"] <- res[[i]]["rd_imp", "length"]
  }
  
  cilength <- colMeans(cilength) # ; round(cilength, 2)
  
  rejrate <- matrix(0, nrow = length(res), ncol = 4) 
  colnames(rejrate) <- c("new", "rtb", "wash", "rd_imp")
  for(i in 1:nrow(rejrate)) {
    
    rejrate[i, "new"] <- res[[i]]["new", "rejected"] 
    rejrate[i, "rtb"] <- res[[i]]["rtb", "rejected"] 
    rejrate[i, "wash"] <- res[[i]]["wash", "rejected"] 
    rejrate[i, "rd_imp"] <- res[[i]]["rd_imp", "rejected"]
  }
  
  rejrate <- colMeans(rejrate) # ; round(rejrate, 2)
  
  output <- data.frame(BETA_HYP = round(beta_hy, 3), 
                       BETA_TP = round(beta_tp, 3), 
                       BIAS = round(bias, 3), 
                       RMSE = round(rmse, 3),  
                       REJRATE = round(rejrate, 3), 
                       COVRATE = round(covrate, 2), 
                       CI95LENGTH = round(cilength, 3),
                       row.names = c("NEW", "RTB", "WASHOUT", "RDIMPUTE"))
  
  output[c("RTB", "WASHOUT", "RDIMPUTE", "NEW"), ]
}

##############################################################################
###
##############################################################################

# n = 200 and pi = 0.5
outtab(res04, beta_hy04, beta_tp04) ; outtab(res05, beta_hy05, beta_tp05) ; outtab(res06, beta_hy06, beta_tp06)
outtab(res01, beta_hy01, beta_tp01) ; outtab(res02, beta_hy02, beta_tp02) ; outtab(res03, beta_hy03, beta_tp03)

# n = 100 and pi = 0.5
outtab(res10, beta_hy10, beta_tp10) ; outtab(res11, beta_hy11, beta_tp11) ; outtab(res12, beta_hy12, beta_tp12)
outtab(res07, beta_hy07, beta_tp07) ; outtab(res08, beta_hy08, beta_tp08) ; outtab(res09, beta_hy09, beta_tp09)

# n = 200 and pi = 0.3
outtab(res16, beta_hy16, beta_tp16) ; outtab(res17, beta_hy17, beta_tp17) ; outtab(res18, beta_hy18, beta_tp18)
outtab(res13, beta_hy13, beta_tp13) ; outtab(res14, beta_hy14, beta_tp14) ; outtab(res15, beta_hy15, beta_tp15)

# n = 200 and pi = 0.1
outtab(res22, beta_hy22, beta_tp22) ; outtab(res23, beta_hy23, beta_tp23) ; outtab(res24, beta_hy24, beta_tp24)
outtab(res19, beta_hy19, beta_tp19) ; outtab(res20, beta_hy20, beta_tp20) ; outtab(res21, beta_hy21, beta_tp21)

save.image(file = "simout.RData")