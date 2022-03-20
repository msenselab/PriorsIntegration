## the definition of the function to run baseline model for BR shor and long session
funFitStanBaseline_log <- function(subdat, myBaselineModel){
  #subdat <- sub_exp[[1]]
  library(rstan)
  library(tidyverse)
  library(dplyr)
  library(loo)
  
  NewY_mix_list = NULL
  PredY_mix_list = NULL
  subNo <- unique(subdat$NSub)
  expName <- unique(subdat$Exp)
  
  print(paste0('Start run rstan baseline model on Subject No.',subNo, ' in ', expName))
  
  xmean <- subdat %>% dplyr::group_by(group) %>% dplyr::summarise(targetMean =mean(targetDur), RPMean =mean(RP),  xsd = sd(targetDur))
  data_s<- subdat %>% dplyr::filter(group == 'short')  # short groups 
  data_l <- subdat %>% dplyr::filter(group == 'long')  # long groups 
  data_mix <- subdat %>% dplyr::filter(group == 'mixed')  # mixed groups 
  data_mix$model = 'IP'
  
  n_s=length(data_s$RP)
  n_l=length(data_l$RP)
  n_mix = length(data_mix$RP)
  xsnew <- seq(0.4, 0.8, 0.01)  #41
  xlnew <- seq(1.2, 2.4,0.01)  #121
  xnew <- c(seq(0.4, 0.8,0.01), seq(1.2, 2.4,0.01))   #162
  
  stan_data1 = list(Y_s=data_s$RP, n_s=n_s, X_s = data_s$targetDur,
                    Y_l=data_l$RP, n_l=n_l, X_l = data_l$targetDur, 
                    xmean = xmean$targetMean, 
                    xsd = xmean$xsd,
                    xsnew = xsnew,
                    xlnew = xlnew,
                    n_mix = n_mix, 
                    X_mix = data_mix$targetDur,
                    xnew = xnew
  )  #data passed to stan
  
  
  PredY_s_list <- data_s[c('NSub','targetDur', 'RP','Exp','group')]
  PredY_l_list <- data_l[c('NSub','targetDur', 'RP','Exp','group')]
  
  
  myinits1 <- list(
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1),
    list(mu_p_s_log= 0, mu_p_l_log= 0, sig2_mn =0.1, sig_pr2_s_log =0.1, sig_pr2_l_log =0.1, sig2_t = 0.1)
  )
  
  
  # fit the baseline model
  parameters1 <- c("mu_p_s_log","mu_p_l_log", "sig_pr2_s_log", "sig_pr2_l_log",  "sig2_t", "sig2_mn", 
                   "ynew_s", "ynew_l", "ypred_s",  "ypred_l", "ynew_mix", "ypred_mix")
  
  subfit <- sampling(myBaselineModel, 
                     stan_data1,
                     init=myinits1,
                     iter=12000,
                     chains=8,
                     thin=1,
                     cores = parallel::detectCores())
  
  fitpar <- summary(subfit, pars = parameters1)$summary
  list_of_draws <- rstan::extract(subfit, pars = parameters1)
  log_lik_rlt <- extract_log_lik(subfit)
  loo_1 <- loo(log_lik_rlt)
  waic = waic(log_lik_rlt)
  mu_p_s_log =  mean(list_of_draws$mu_p_s_log)
  mu_p_l_log =  mean(list_of_draws$mu_p_l_log)
  sig_pr2_s_log = mean(list_of_draws$sig_pr2_s_log)
  sig_pr2_l_log = mean(list_of_draws$sig_pr2_l_log)
  sig2_t = mean(list_of_draws$sig2_t)
  sig2_mn = mean(list_of_draws$sig2_mn)
  
  
  #prediction of mix session
  ynew_mix_list <- list_of_draws$ypred_mix
  newy_mix <- matrix(rep(0, 162, 7), nrow = 162, ncol = 7)
  for (i in 1:7){
    for (j in 1:162){
      newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
    }
  } 
  NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
  colnames(NewY_mix_list)  <-c("curDur", "W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
  NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_P_G)
  NewY_mix_list$W_Ds = (1-NewY_mix_list$W_DL_PL) * (1-NewY_mix_list$W_P_G)
  NewY_mix_list$W_P_S = NewY_mix_list$W_L
  NewY_mix_list$W_P_L = NewY_mix_list$W_L
  NewY_mix_list[which(NewY_mix_list$curDur <= 1),"W_P_L"] = 0
  NewY_mix_list[which(NewY_mix_list$curDur >= 1),"W_P_S"] = 0
  NewY_mix_list$NSub =subNo
  NewY_mix_list$exp = expName
  NewY_mix_list$model = 'IP' #independent prior based on short and long session
  NewY_mix_list$group = 'mixed'
  
  pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
  y_pred_mix_list <- list_of_draws$ynew_mix
  for (i in 1:7){
    for (j in 1:n_mix){
      pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
    }
  }
  
  pred_y_mix = data.frame(pred_y_mix)
  colnames(pred_y_mix)  <-c("W_P_G", "mu_r", "sig_r", "predY", "W_PL", "sig2_DL","log_lik")
  pred_y_mix$W_L = pred_y_mix$W_PL * (1-pred_y_mix$W_P_G)
  pred_y_mix$W_Ds = (1- pred_y_mix$W_PL) * (1-pred_y_mix$W_P_G)
  pred_y_mix$W_P_S = pred_y_mix$W_L
  pred_y_mix$W_P_L = pred_y_mix$W_L
  pred_y_mix[which(pred_y_mix$curDur <= 1),"W_P_L"] = 0
  pred_y_mix[which(pred_y_mix$curDur >= 1),"W_P_S"] = 0
  PredY_mix_list <- data_mix[c('NSub','targetDur', 'RP','Exp','group', 'model')]
  PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
  
  #prediction of short session
  newy_s <- matrix(rep(0, 41, 4), nrow = 41, ncol = 4)
  ynew_s_list <- list_of_draws$ypred_s
  
  for (i in 1:4){
    for (j in 1:41){
      newy_s[j,i] <-  mean(ynew_s_list[,j,i] )
    }
  }
  NewY_s_list <- data.frame(cbind(xsnew, newy_s))
  colnames(NewY_s_list)  <-c("curDur", "wp", "mu_r", "sig_r", "predY")
  NewY_s_list$NSub =subNo
  NewY_s_list$exp = expName
  NewY_s_list$model = 'baseline'
  NewY_s_list$group = 'short'
  
  #prediction of long session
  newy_l <- matrix(rep(0, 121, 4), nrow = 121, ncol = 4)
  ynew_l_list <- list_of_draws$ypred_l
  for (i in 1:4){
    for (j in 1:121){
      newy_l[j,i] <-  mean(ynew_l_list[,j,i] )
    }
  }
  NewY_l_list <- data.frame(cbind(xlnew, newy_l))
  colnames(NewY_l_list)  <-c("curDur", "wp", "mu_r", "sig_r", "predY")
  NewY_l_list$NSub =subNo
  NewY_l_list$exp = expName
  NewY_l_list$model = 'baseline'
  NewY_l_list$group = 'long'
  
  pred_y_s <- matrix(rep(0, n_s, 4), nrow = n_s, ncol = 4)
  y_pred_s_list <- list_of_draws$ynew_s
  #print(dim(predRP_list))
  for (i in 1:4){
    for (j in 1:n_s){
      pred_y_s[j,i] <-  mean(y_pred_s_list[,j,i] )
    }
  }
  colnames(pred_y_s)  <-c("wp", "mu_r", "sig_r", "predY")
  PredY_s_list <- cbind(PredY_s_list, pred_y_s)
  
  
  pred_y_l <- matrix(rep(0, n_l, 4), nrow = n_l, ncol = 4)
  y_pred_l_list <- list_of_draws$ynew_l
  for (i in 1:4){
    for (j in 1:n_l){
      pred_y_l[j,i] <-  mean(y_pred_l_list[,j,i] )
    }
  }
  colnames(pred_y_l)  <-c("wp", "mu_r", "sig_r", "predY")
  PredY_l_list <- cbind(PredY_l_list, pred_y_l)
  
  Baypar = data.frame(
    Nsub = subNo,
    Exp = expName,
    model = 'Baseline',
    mu_p_s_log = mu_p_s_log,
    mu_p_l_log = mu_p_l_log, 
    sig_pr2_s_log = sig_pr2_s_log,
    sig_pr2_l_log = sig_pr2_l_log,
    sig2_mn = sig2_mn,
    sig2_t = sig2_t,
    looic = loo_1$looic,
    p_loo = loo_1$p_loo,
    elpd_loo = loo_1$elpd_loo,
    se_looic = loo_1$se_looic,
    se_p_loo = loo_1$se_p_loo,
    waic = waic$waic,
    p_waic =waic$p_waic,
    se_waic = waic$se_waic,
    se_p_waic = waic$se_p_waic,
    elpd_waic = waic$elpd_waic
  )
  
  return(list("Baypar" = Baypar, "NewY_s_list" = NewY_s_list,
              "NewY_l_list" = NewY_l_list, 
              "PredY_s_list" = PredY_s_list, "PredY_l_list" = PredY_l_list, 
              "PredY_mix_list" = PredY_mix_list, "NewY_mix_list" = NewY_mix_list))
}


funFitStanBaseline_linear <- function(subdat, myBaselineModel){
  #subdat <- sub_exp[[1]]
  library(rstan)
  library(tidyverse)
  library(dplyr)
  library(loo)
  
  NewY_mix_list = NULL
  PredY_mix_list = NULL
  subNo <- unique(subdat$NSub)
  expName <- unique(subdat$Exp)
  
  print(paste0('Start run rstan baseline model on Subject No.',subNo, ' in ', expName))
  xmean <- subdat %>% dplyr::group_by(group) %>% dplyr::summarise(targetMean =mean(targetDur), RPMean =mean(RP),  xsd = sd(targetDur))
  data_s<- subdat %>% dplyr::filter(group == 'short')  # short session 
  data_l <- subdat %>% dplyr::filter(group == 'long')  # long session 
  data_mix <- subdat %>% dplyr::filter(group == 'mixed')  # mixed session 
  data_mix$model = 'IP'
  
  n_s=length(data_s$RP)
  n_l=length(data_l$RP)
  n_mix = length(data_mix$RP)
  xsnew <- seq(0.4, 0.8, 0.01)  #41
  xlnew <- seq(1.2, 2.4,0.01)  #121
  xnew <- c(seq(0.4, 0.8,0.01), seq(1.2, 2.4,0.01))   #162
  
  stan_data1 = list(Y_s=data_s$RP, n_s=n_s, X_s = data_s$targetDur,
                    Y_l=data_l$RP, n_l=n_l, X_l = data_l$targetDur, 
                    xmean = xmean$targetMean, 
                    xsd = xmean$xsd,
                    xsnew = xsnew,
                    xlnew = xlnew,
                    n_mix = n_mix, 
                    X_mix = data_mix$targetDur,
                    xnew = xnew
  )  #data passed to stan
  
  
  PredY_s_list <- data_s[c('NSub','targetDur', 'RP','Exp','group')]
  PredY_l_list <- data_l[c('NSub','targetDur', 'RP','Exp','group')]
  
  m_sdur = mean(data_s$targetDur)
  m_ldur = mean(data_l$targetDur)
  myinits1 <- list(
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1),
    list(mu_p_s= m_sdur, mu_p_l= m_ldur, sig2_mn =0.1, sig_pr2_s =0.1, sig_pr2_l =0.1, sig_t = 0.1)
  )
  
  
  # fit the baseline model
  parameters1 <- c("mu_p_s","mu_p_l", "sig_pr2_s", "sig_pr2_l",  "sig_t", "sig2_mn", 
                   "ynew_s", "ynew_l", "ypred_s",  "ypred_l", "ynew_mix", "ypred_mix")
  
  subfit <- sampling(myBaselineModel, 
                     stan_data1,
                     init=myinits1,
                     iter=12000,
                     chains=8,
                     thin=1,
                     cores = parallel::detectCores())
  
  fitpar <- summary(subfit, pars = parameters1)$summary
  list_of_draws <- rstan::extract(subfit, pars = parameters1)
  log_lik_rlt <- extract_log_lik(subfit)
  loo_1 <- loo(log_lik_rlt)
  waic = waic(log_lik_rlt)
  mu_p_s =  mean(list_of_draws$mu_p_s)
  mu_p_l =  mean(list_of_draws$mu_p_l)
  sig_pr2_s = mean(list_of_draws$sig_pr2_s)
  sig_pr2_l = mean(list_of_draws$sig_pr2_l)
  sig_t = mean(list_of_draws$sig_t)
  sig2_mn = mean(list_of_draws$sig2_mn)
  
  
  #prediction of mix session
  ynew_mix_list <- list_of_draws$ypred_mix
  newy_mix <- matrix(rep(0, 162, 7), nrow = 162, ncol = 7)
  for (i in 1:7){
    for (j in 1:162){
      newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
    }
  } 
  NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
  colnames(NewY_mix_list)  <-c("curDur", "W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
  NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_P_G)
  NewY_mix_list$W_Ds = (1-NewY_mix_list$W_DL_PL) * (1-NewY_mix_list$W_P_G)
  NewY_mix_list$W_P_S = NewY_mix_list$W_L
  NewY_mix_list$W_P_L = NewY_mix_list$W_L
  NewY_mix_list[which(NewY_mix_list$curDur <= 1),"W_P_L"] = 0
  NewY_mix_list[which(NewY_mix_list$curDur >= 1),"W_P_S"] = 0
  NewY_mix_list$NSub =subNo
  NewY_mix_list$exp = expName
  NewY_mix_list$model = 'IP' #independent prior based on short and long session
  NewY_mix_list$group = 'mixed'
  
  pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
  y_pred_mix_list <- list_of_draws$ynew_mix
  for (i in 1:7){
    for (j in 1:n_mix){
      pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
    }
  }
  
  pred_y_mix = data.frame(pred_y_mix)
  colnames(pred_y_mix)  <-c("W_P_G", "mu_r", "sig_r", "predY", "W_PL", "sig2_DL","log_lik")
  pred_y_mix$W_L = pred_y_mix$W_PL * (1-pred_y_mix$W_P_G)
  pred_y_mix$W_Ds = (1- pred_y_mix$W_PL) * (1-pred_y_mix$W_P_G)
  pred_y_mix$W_P_S = pred_y_mix$W_L
  pred_y_mix$W_P_L = pred_y_mix$W_L
  pred_y_mix[which(pred_y_mix$curDur <= 1),"W_P_L"] = 0
  pred_y_mix[which(pred_y_mix$curDur >= 1),"W_P_S"] = 0
  PredY_mix_list <- data_mix[c('NSub','targetDur', 'RP','Exp','group', 'model')]
  PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
  
  
  #prediction of short session
  newy_s <- matrix(rep(0, 41, 4), nrow = 41, ncol = 4)
  ynew_s_list <- list_of_draws$ypred_s
  
  for (i in 1:4){
    for (j in 1:41){
      newy_s[j,i] <-  mean(ynew_s_list[,j,i] )
    }
  }
  NewY_s_list <- data.frame(cbind(xsnew, newy_s))
  colnames(NewY_s_list)  <-c("curDur", "wp", "mu_r", "sig_r", "predY")
  NewY_s_list$NSub =subNo
  NewY_s_list$exp = expName
  NewY_s_list$model = 'baseline'
  NewY_s_list$group = 'short'
  
  #prediction of long session
  newy_l <- matrix(rep(0, 121, 4), nrow = 121, ncol = 4)
  ynew_l_list <- list_of_draws$ypred_l
  for (i in 1:4){
    for (j in 1:121){
      newy_l[j,i] <-  mean(ynew_l_list[,j,i] )
    }
  }
  NewY_l_list <- data.frame(cbind(xlnew, newy_l))
  colnames(NewY_l_list)  <-c("curDur", "wp", "mu_r", "sig_r", "predY")
  NewY_l_list$NSub =subNo
  NewY_l_list$exp = expName
  NewY_l_list$model = 'baseline'
  NewY_l_list$group = 'long'
  
  pred_y_s <- matrix(rep(0, n_s, 4), nrow = n_s, ncol = 4)
  y_pred_s_list <- list_of_draws$ynew_s
  #print(dim(predRP_list))
  for (i in 1:4){
    for (j in 1:n_s){
      pred_y_s[j,i] <-  mean(y_pred_s_list[,j,i] )
    }
  }
  colnames(pred_y_s)  <-c("wp", "mu_r", "sig_r", "predY")
  PredY_s_list <- cbind(PredY_s_list, pred_y_s)
  
  
  pred_y_l <- matrix(rep(0, n_l, 4), nrow = n_l, ncol = 4)
  y_pred_l_list <- list_of_draws$ynew_l
  for (i in 1:4){
    for (j in 1:n_l){
      pred_y_l[j,i] <-  mean(y_pred_l_list[,j,i] )
    }
  }
  colnames(pred_y_l)  <-c("wp", "mu_r", "sig_r", "predY")
  PredY_l_list <- cbind(PredY_l_list, pred_y_l)
  
  Baypar = data.frame(
    Nsub = subNo,
    Exp = expName,
    model = 'Baseline',
    mu_p_s = mu_p_s,
    mu_p_l = mu_p_l, 
    sig_pr2_s = sig_pr2_s,
    sig_pr2_l = sig_pr2_l,
    sig2_mn = sig2_mn,
    sig_t = sig_t,
    looic = loo_1$looic,
    p_loo = loo_1$p_loo,
    elpd_loo = loo_1$elpd_loo,
    se_looic = loo_1$se_looic,
    se_p_loo = loo_1$se_p_loo,
    waic = waic$waic,
    p_waic =waic$p_waic,
    se_waic = waic$se_waic,
    se_p_waic = waic$se_p_waic,
    elpd_waic = waic$elpd_waic
  )
  
  PredY_mix_list$part= partNo
  return(list("Baypar" = Baypar, "NewY_s_list" = NewY_s_list,
              "NewY_l_list" = NewY_l_list, 
              "PredY_s_list" = PredY_s_list, "PredY_l_list" = PredY_l_list, 
              "PredY_mix_list" = PredY_mix_list, "NewY_mix_list" = NewY_mix_list))
}


# funFitStanSub_linear <- function(data_mix, myrstanModel, partNo){
#   #data_mix = subdata_mix%>%filter(id < total/2)
#   library(rstan)
#   library(tidyverse)
#   library(dplyr)
#   library(loo)
#   
#   modelname <-  unique(data_mix$model)
#   subNo <- unique(data_mix$NSub)
#   expName <- unique(data_mix$Exp)
#   AllDat_Bayparlist_BR<- read.csv("RSTANMODELS/models/Baseline_linear/AllDat_Bayparlist_BR.csv")
#   
#   BR_par = AllDat_Bayparlist_BR%>% filter(Nsub == subNo, Exp ==expName)
#   
#   print(paste0('Start run rstan model ', modelname,' on mixed session of Subject No.',subNo, ' in ', expName, ':', partNo ))
#   
#   durList <- sort(unique(data_mix$targetDur))
#   for(i in 1: length(durList)){
#     data_mix[which(data_mix$targetDur == durList[i]),"DurIdx"] = i
#   }
#   
#   n_mix = length(data_mix$RP)
#   xnew <- c(seq(0.4, 0.8,0.01), seq(1.2, 2.4,0.01))   #162
#   mu_p_s = BR_par$mu_p_s
#   mu_p_l = BR_par$mu_p_l
#   sig_pr2_s = BR_par$sig_pr2_s
#   sig_pr2_l = BR_par$sig_pr2_l
#   
#   stan_data2 = list( mu_p_s = mu_p_s,
#                      mu_p_l = mu_p_l, 
#                      sig_pr2_s = sig_pr2_s,
#                      sig_pr2_l = sig_pr2_l,
#                      sig2_mn = BR_par$sig2_mn,
#                      sig_t = BR_par$sig_t,
#                      Y_mix = data_mix$RP, 
#                      n_mix = n_mix, 
#                      X_mix = data_mix$targetDur,
#                      X_idx = data_mix$DurIdx,
#                      xnew = xnew
#   )  #data passed to stan
#   
#   parameters2 <- c( "mu_p_g", "sig2_p_g", "ynew_mix", "ypred_mix")
#   myinits2 <- list(
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1)
#   )
#   
#   
#   
#   if (modelname == 'BPM'){
#     parameters2 <- c("mp", "ynew_mix") 
#     myinits2 <- list(list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5))
#   }
#   
#   if (modelname  == 'DIM'){
#     parameters2 <- c( "mu_p_g", "sig2_p_g", "phi", "ynew_mix", "ypred_mix") 
#   }
#   
#   if (modelname == 'IP') {
#     parameters2 <- c("mu_p_s", "mu_p_l","sig_pr2_s", "sig_pr2_l", "ynew_mix", "ypred_mix") 
#     myinits2 <- list(
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),    
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1)
#     )
#   }
#   
#   
#   subfit2 <- sampling(myrstanModel, 
#                       stan_data2,
#                       init=myinits2,
#                       iter=4000,
#                       chains=8,
#                       thin=1,
#                       control = list(adapt_delta = 0.99,
#                                      max_treedepth = 15),
#                       cores = parallel::detectCores())
#   
#   
#   fitpar <- summary(subfit2, pars = parameters2)$summary
#   list_of_draws <- rstan::extract(subfit2, pars = parameters2)
#   phi = 0
#   mp = NULL
#   phi_list = NULL
#   Baypar = NULL
#   NewY_mix_list = NULL
#   PredY_mix_list = NULL
#   
#   log_lik_rlt <- extract_log_lik(subfit2)
#   loo_1 <- loo(log_lik_rlt)
#   waic = waic(log_lik_rlt)
#   
#   if(modelname == 'IP'){
#     mu_p_s_BR =  mean(list_of_draws$mu_p_s)
#     mu_p_l_BR =  mean(list_of_draws$mu_p_l)
#     sig_pr2_s_BR = mean(list_of_draws$sig_pr2_s)
#     sig_pr2_l_BR = mean(list_of_draws$sig_pr2_l)
#   }else{
#     mu_p_s_BR =  mu_p_s
#     mu_p_l_BR =  mu_p_l
#     sig_pr2_s_BR = sig_pr2_s
#     sig_pr2_l_BR = sig_pr2_l
#   }
#   
#   if(modelname == 'BPM'){
#     mp0 =  mean(list_of_draws$mp)
#     
#     mp <-  c(mp0, 1-mp0, 0, subNo, expName, modelname)
#     mp = data.frame(mp)
#     
#     pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
#     y_pred_mix_list <- list_of_draws$ypred_mix
#     for (i in 1:7){
#       for (j in 1:n_mix){
#         pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#       }
#     }
#     pred_y_mix = data.frame(pred_y_mix)
#     colnames(pred_y_mix)  <-c("predY_s", "predY_l", "predY", "mu_r", "wp_ps", "wp_pl", "log_lik")
#     
#     PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model', 'DurIdx')]
#     PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
#     
#     PredY_mix_list$mp = mp0
#     PredY_mix_list$W_P_S = mp0*PredY_mix_list$wp_ps
#     PredY_mix_list$W_P_L = (1-mp0)* PredY_mix_list$wp_pl
#     PredY_mix_list$W_Ds = mp0*(1-PredY_mix_list$wp_ps)+ (1-mp0)*(1-PredY_mix_list$wp_pl)
#     PredY_mix_list$W_P_G = 0
#     PredY_mix_list$group ='mixed'
#     NewY_mix_list$part= partNo
#     PredY_mix_list$sig_r =0. # Question: how to calculate sig_r in BPM model?
#     
#     mu_p_g = 1
#     sig2_p_g = 0.1
#   }else{
#     if(modelname == 'IP'){
#       mu_p_g = 1
#       sig2_p_g = 0.1
#     }else{
#       mu_p_g =  mean(list_of_draws$mu_p_g)
#       sig2_p_g = mean(list_of_draws$sig2_p_g)
#     }
#     
#     if(modelname == 'DIM'){
#       phi = mean(list_of_draws$phi)
#       ynew_mix_list <- list_of_draws$ynew_mix
#       newy_mix <- matrix(rep(0, 162, 9), nrow = 162, ncol = 9)
#       for (i in 1:9){
#         for (j in 1:162){
#           newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
#         }
#       }
#       
#       NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#       # 1st col //weight of D_G in De  "W_DG"
#       # 2nd col mu_r # 3th col sig_r,  #4th col predY
#       # 5th col weight of local prior PL in DL  "W_DL_PL"  
#       # 6th col weight of global prior in DG  "W_DG_PG"
#       # 7th col variance of DL   "sig2_DL" 
#       # 8th col variance of DG  "sig2_DG"
#       # 9th col likelihood
#       colnames(NewY_mix_list)  <-c("curDur","W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
#       #weight of global prior in De
#       NewY_mix_list$W_P_G =  NewY_mix_list$W_DG_PG * NewY_mix_list$W_DG
#       #weight of local prior in De
#       NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_DG)
#       #weight of Ds in De 
#       NewY_mix_list$W_Ds =(1-NewY_mix_list$W_DL_PL) *(1-NewY_mix_list$W_DG) + NewY_mix_list$W_DG* (1- NewY_mix_list$W_DG_PG)
#       
#     }
#     
#     if(modelname %in% c('LGM', 'PIM', 'IP')){
#       ynew_mix_list <- list_of_draws$ynew_mix
#       newy_mix <- matrix(rep(0, 162, 7), nrow = 162, ncol = 7)
#       for (i in 1:7){
#         for (j in 1:162){
#           newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
#         }
#       }
#       
#       if(modelname%in% c('LGM','IP')){
#         # 1st col weight of global prior in De  "W_P_G"
#         # 2nd col mu_r # 3th col sig_r,  #4th col predY
#         # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
#         # 6th col variance of DL   "sig2_DL" 
#         # 7th col likelihood
#         NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#         colnames(NewY_mix_list)  <-c("curDur", "W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
#         NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_P_G)
#         NewY_mix_list$W_Ds = (1-NewY_mix_list$W_DL_PL) * (1-NewY_mix_list$W_P_G)
#         
#       }else if(modelname == 'PIM'){
#         NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#         # 1st col weight of Ds in De  "W_Ds"
#         # 2nd col mu_r # 3th col sig_r,  #4th col predY
#         # 5th col weight of local prior PL in integrated prior PI  "W_PI_PL" 
#         # 6th col variance of PI   "sig2_PI" 
#         # 7th col likelihood
#         colnames(NewY_mix_list)  <-c("curDur","W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
#         # weight of global prior in De
#         NewY_mix_list$W_P_G = (1-NewY_mix_list$W_Ds) * (1-NewY_mix_list$W_PI_PL)
#         # weight of local prior in De
#         NewY_mix_list$W_L = (1- NewY_mix_list$W_Ds) * NewY_mix_list$W_PI_PL
#       }
#       
#       if(modelname %in% c('DIM', 'LGM', 'PIM', 'IP')){
#         NewY_mix_list$W_P_S = NewY_mix_list$W_L
#         NewY_mix_list$W_P_L = NewY_mix_list$W_L
#         NewY_mix_list[which(NewY_mix_list$curDur <= 1),"W_P_L"] = 0
#         NewY_mix_list[which(NewY_mix_list$curDur >= 1),"W_P_S"] = 0
#       }
#       
#       
#       NewY_mix_list$NSub =subNo
#       NewY_mix_list$exp = expName
#       NewY_mix_list$model = modelname
#       NewY_mix_list$group = 'mixed'
#       NewY_mix_list$part= partNo
#       
#       pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
#       y_pred_mix_list <- list_of_draws$ypred_mix
#       for (i in 1:7){
#         for (j in 1:n_mix){
#           pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#         }
#       }
#       
#       pred_y_mix = data.frame(pred_y_mix)
#       if(modelname %in% c('LGM', 'IP')){
#         # 1st col weight of global prior in De  "W_P_G"
#         # 2nd col mu_r # 3th col sig_r,  #4th col predY
#         # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
#         # 6th col variance of DL   "sig2_DL" 
#         # 7th col likelihood
#         colnames(pred_y_mix)  <- c("W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL", "sig2_DL","log_lik")
#         pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_P_G)
#         pred_y_mix$W_Ds = (1- pred_y_mix$W_DL_PL) * (1-pred_y_mix$W_P_G)
#       }else if(modelname == 'PIM'){
#         colnames(pred_y_mix)  <-c("W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
#         # weight of global prior in De
#         pred_y_mix$W_P_G = (1-pred_y_mix$W_Ds) * (1-pred_y_mix$W_PI_PL)
#         # weight of local prior in De
#         pred_y_mix$W_L = (1- pred_y_mix$W_Ds) * pred_y_mix$W_PI_PL
#       }
#     }
#     
#     
#     #for DIM results
#     if(modelname %in% c('DIM')){
#       pred_y_mix <- matrix(rep(0, n_mix, 9), nrow = n_mix, ncol = 9)
#       y_pred_mix_list <- list_of_draws$ypred_mix
#       for (i in 1:9){
#         for (j in 1:n_mix){
#           pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#         }
#       }
#       pred_y_mix = data.frame(pred_y_mix)
#       # 1st col weight of D_G in De  "W_DG"
#       # 2nd col mu_r # 3th col sig_r,  #4th col predY
#       # 5th col weight of local prior PL in DL  "W_DL_PL"  
#       # 6th col weight of global prior in DG  "W_DG_PG"
#       # 7th col variance of DL   "sig2_DL" 
#       # 8th col variance of DG  "sig2_DG"
#       # 9th col likelihood
#       colnames(pred_y_mix)  <-c("W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
#       #weight of global prior in De
#       pred_y_mix$W_P_G =  pred_y_mix$W_DG_PG * pred_y_mix$W_DG
#       #weight of local prior in De
#       pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_DG)
#       #weight of Ds in De 
#       pred_y_mix$W_Ds = (1-pred_y_mix$W_DL_PL) *(1-pred_y_mix$W_DG) + pred_y_mix$W_DG* (1-pred_y_mix$W_DG_PG)
#     } 
#     
#     if(modelname %in% c('DIM', 'LGM', 'PIM', 'IP')){
#       pred_y_mix$W_P_S = pred_y_mix$W_L
#       pred_y_mix$W_P_L = pred_y_mix$W_L
#       pred_y_mix[which(pred_y_mix$curDur <= 1),"W_P_L"] = 0
#       pred_y_mix[which(pred_y_mix$curDur >= 1),"W_P_S"] = 0
#     }
#     PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model')]
#     PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
#   }
#   
#   Baypar = data.frame(
#     Nsub = subNo,
#     Exp = expName,
#     model = modelname,
#     phi = phi,
#     mu_p_s = mu_p_s,
#     mu_p_l = mu_p_l, 
#     mu_p_g = mu_p_g,
#     sig_pr2_s = sig_pr2_s,
#     sig_pr2_l = sig_pr2_l,
#     sig2_p_g = sig2_p_g,
#     sig2_mn = BR_par$sig2_mn,
#     sig_t = BR_par$sig_t,
#     looic = loo_1$looic,
#     p_loo = loo_1$p_loo,
#     elpd_loo = loo_1$elpd_loo,
#     se_looic = loo_1$se_looic,
#     se_p_loo = loo_1$se_p_loo,
#     waic = waic$waic,
#     p_waic =waic$p_waic,
#     se_waic = waic$se_waic,
#     se_p_waic = waic$se_p_waic,
#     elpd_waic = waic$elpd_waic
#   )
#   
#   Baypar_BR = data.frame(
#     Nsub = subNo,
#     Exp = expName,
#     model = modelname,
#     mu_p_s = mu_p_s,
#     mu_p_l = mu_p_l, 
#     sig_pr2_s = sig_pr2_s,
#     sig_pr2_l = sig_pr2_l,
#     sig2_mn = BR_par$sig2_mn,
#     sig_t = BR_par$sig_t,
#     mu_p_s_BR =  mu_p_s_BR,
#     mu_p_l_BR =  mu_p_l_BR,
#     sig_pr2_s_BR = sig_pr2_s_BR,
#     sig_pr2_l_BR = sig_pr2_l_BR
#   )
#   PredY_mix_list$part = partNo
#   return(list("Baypar" = Baypar, "NewY_mix_list" = NewY_mix_list, "PredY_mix_list" = PredY_mix_list, 
#               "mp" = mp, "phi_list" = phi_list, "Baypar_BR" = Baypar_BR))
# }
# 
# 
# ## The definition of the function to run Rstan model with sub data
# funFitStanSub_log <- function(data_mix, myrstanModel, partNo){
#   modelname <-  unique(data_mix$model)
#   subNo <- unique(data_mix$NSub)
#   expName <- unique(data_mix$Exp)
#   AllDat_Bayparlist_BR<- read.csv("RSTANMODELS/models/Baseline_log/AllDat_Bayparlist_BR.csv")
#   BR_par = AllDat_Bayparlist_BR%>% filter(Nsub == subNo, Exp ==expName)
#   print(paste0('Start run rstan model ', modelname,' on mixed session of Subject No.',subNo, ' in ', expName, ':', partNo ))
#   
#   durList <- sort(unique(data_mix$targetDur))
#   for(i in 1: length(durList)){
#     data_mix[which(data_mix$targetDur == durList[i]),"DurIdx"] = i
#   }
#   
#   n_mix = length(data_mix$RP)
#   xsnew <- seq(0.4, 0.8, 0.01)  #41
#   xlnew <- seq(1.2, 2.4,0.01)  #121
#   xnew <- c(seq(0.4, 0.8,0.01), seq(1.2, 2.4,0.01))   #162
#   
#   mu_p_s_log = BR_par$mu_p_s_log
#   mu_p_l_log = BR_par$mu_p_l_log
#   sig_pr2_s_log = BR_par$sig_pr2_s_log
#   sig_pr2_l_log = BR_par$sig_pr2_l_log
#   
#   stan_data2 = list( mu_p_s_log = mu_p_s_log,
#                      mu_p_l_log = mu_p_l_log, 
#                      sig_pr2_s_log = sig_pr2_s_log,
#                      sig_pr2_l_log = sig_pr2_l_log,
#                      sig2_mn = BR_par$sig2_mn,
#                      sig2_t = BR_par$sig2_t,
#                      Y_mix = data_mix$RP, 
#                      n_mix = n_mix, 
#                      X_mix = data_mix$targetDur,
#                      X_idx = data_mix$DurIdx,
#                      xnew = xnew)  #data passed to stan
#   
#   parameters2 <- c("mu_p_g", "sig2_p_g", "ynew_mix", "ypred_mix")
#   myinits2 <- list(
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),    
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1),
#     list(mu_p_g=1, sig2_p_g = 0.1)
#   )
#   
#   
#   if (modelname == 'BPM'){
#     parameters2 <- c("mp", "ynew_mix") 
#     myinits2 <- list(list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5),
#                      list(mp = 0.5))
#   }
#   
#   if (modelname  == 'DIM'){
#     parameters2 <- c( "mu_p_g", "sig2_p_g", "phi", "ynew_mix", "ypred_mix") 
#   }
#   
#   if (modelname == 'IP') {
#     parameters2 <- c("mu_p_s", "mu_p_l","sig_pr2_s", "sig_pr2_l", "ynew_mix", "ypred_mix") 
#     myinits2 <- list(
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),    
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
#       list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1)
#     )
#   }
#   
#   
#   subfit2 <- sampling(myrstanModel, 
#                       stan_data2,
#                       init=myinits2,
#                       iter=4000,
#                       chains=8,
#                       thin=1,
#                       control = list(adapt_delta = 0.9,
#                                      max_treedepth = 15),
#                       cores = parallel::detectCores())
#   
#   fitpar <- summary(subfit2, pars = parameters2)$summary
#   list_of_draws <- rstan::extract(subfit2, pars = parameters2)
#   phi = 0
#   mp = NULL
#   phi_list = NULL
#   Baypar = NULL
#   NewY_mix_list = NULL
#   PredY_mix_list = NULL
#   
#   log_lik_rlt <- extract_log_lik(subfit2)
#   loo_1 <- loo(log_lik_rlt)
#   waic = waic(log_lik_rlt)
#   mu_p_g =  0
#   sig2_p_g = 0   
#   if(modelname == 'IP'){
#     mu_p_s_log_BR =  mu_p_s_log
#     mu_p_l_log_BR =  mu_p_l_log
#     sig_pr2_s_log_BR = sig_pr2_s_log
#     sig_pr2_l_log_BR = sig_pr2_l_log
#     mu_p_s_log_IR =  mean(list_of_draws$mu_p_s_log)
#     mu_p_l_log_IR =  mean(list_of_draws$mu_p_l_log)
#     sig_pr2_s_log_IR = mean(list_of_draws$sig_pr2_s_log)
#     sig_pr2_l_log_IR = mean(list_of_draws$sig_pr2_l_log)
#   }else{
#     mu_p_s_log_BR =  mu_p_s_log
#     mu_p_l_log_BR =  mu_p_l_log
#     sig_pr2_s_log_BR = sig_pr2_s_log
#     sig_pr2_l_log_BR = sig_pr2_l_log
#     mu_p_s_log_IR =  0
#     mu_p_l_log_IR =  0
#     sig_pr2_s_log_IR = 0
#     sig_pr2_l_log_IR = 0
#     if(modelname != 'BPM'){
#       mu_p_g =  mean(list_of_draws$mu_p_g)
#       sig2_p_g = mean(list_of_draws$sig2_p_g)
#     }
#   }
#   
#   
#   if(modelname == 'BPM'){
#     mp0 =  mean(list_of_draws$mp)
#     mp <-  c(mp0, 1-mp0, 0, subNo, expName, modelname)
#     mp = data.frame(mp)
#     
#     pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
#     y_pred_mix_list <- list_of_draws$ypred_mix
#     for (i in 1:7){
#       for (j in 1:n_mix){
#         pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#       }
#     }
#     pred_y_mix = data.frame(pred_y_mix)
#     colnames(pred_y_mix)  <-c("predY_s", "predY_l", "predY", "mu_r", "wp_ps", "wp_pl", "log_lik")
#     
#     PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model', 'DurIdx')]
#     PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
#     
#     PredY_mix_list$mp = mp0
#     PredY_mix_list$W_P_S = mp0*PredY_mix_list$wp_ps
#     PredY_mix_list$W_P_L = (1-mp0)* PredY_mix_list$wp_pl
#     PredY_mix_list$W_Ds = mp0*(1-PredY_mix_list$wp_ps)+ (1-mp0)*(1-PredY_mix_list$wp_pl)
#     PredY_mix_list$W_P_G = 0
#     PredY_mix_list$group ='mixed'
#     NewY_mix_list$part= partNo
#     PredY_mix_list$sig_r =0. # Question: how to calculate sig_r in BPM model?
#   }
#   
#   if(modelname == 'DIM'){
#     phi = mean(list_of_draws$phi)
#     ynew_mix_list <- list_of_draws$ynew_mix
#     newy_mix <- matrix(rep(0, 162, 9), nrow = 162, ncol = 9)
#     for (i in 1:9){
#       for (j in 1:162){
#         newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
#       }
#     }
#     
#     NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#     # 1st col //weight of D_G in De  "W_DG"
#     # 2nd col mu_r # 3th col sig_r,  #4th col predY
#     # 5th col weight of local prior PL in DL  "W_DL_PL"  
#     # 6th col weight of global prior in DG  "W_DG_PG"
#     # 7th col variance of DL   "sig2_DL" 
#     # 8th col variance of DG  "sig2_DG"
#     # 9th col likelihood
#     colnames(NewY_mix_list)  <-c("curDur","W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
#     #weight of global prior in De
#     NewY_mix_list$W_P_G =  NewY_mix_list$W_DG_PG * NewY_mix_list$W_DG
#     #weight of local prior in De
#     NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_DG)
#     #weight of Ds in De 
#     NewY_mix_list$W_Ds =(1-NewY_mix_list$W_DL_PL) *(1-NewY_mix_list$W_DG) + NewY_mix_list$W_DG* (1- NewY_mix_list$W_DG_PG)
#     
#     pred_y_mix <- matrix(rep(0, n_mix, 9), nrow = n_mix, ncol = 9)
#     y_pred_mix_list <- list_of_draws$ypred_mix
#     for (i in 1:9){
#       for (j in 1:n_mix){
#         pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#       }
#     }
#     pred_y_mix = data.frame(pred_y_mix)
#     # 1st col weight of D_G in De  "W_DG"
#     # 2nd col mu_r # 3th col sig_r,  #4th col predY
#     # 5th col weight of local prior PL in DL  "W_DL_PL"  
#     # 6th col weight of global prior in DG  "W_DG_PG"
#     # 7th col variance of DL   "sig2_DL" 
#     # 8th col variance of DG  "sig2_DG"
#     # 9th col likelihood
#     colnames(pred_y_mix)  <-c("W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
#     #weight of global prior in De
#     pred_y_mix$W_P_G =  pred_y_mix$W_DG_PG * pred_y_mix$W_DG
#     #weight of local prior in De
#     pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_DG)
#     #weight of Ds in De 
#     pred_y_mix$W_Ds = (1-pred_y_mix$W_DL_PL) *(1-pred_y_mix$W_DG) + pred_y_mix$W_DG* (1-pred_y_mix$W_DG_PG)
#     
#   }
#   
#   if(modelname %in% c('LGM', 'PIM', 'IP')){
#     ynew_mix_list <- list_of_draws$ynew_mix
#     newy_mix <- matrix(rep(0, 162, 7), nrow = 162, ncol = 7)
#     for (i in 1:7){
#       for (j in 1:162){
#         newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
#       }
#     }
#     
#     pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
#     y_pred_mix_list <- list_of_draws$ypred_mix
#     for (i in 1:7){
#       for (j in 1:n_mix){
#         pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
#       }
#     }
#     pred_y_mix = data.frame(pred_y_mix)
#     
#     if(modelname%in% c('LGM','IP')){
#       NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#       # 1st col weight of global prior in De  "W_P_G"
#       # 2nd col mu_r # 3th col sig_r,  #4th col predY
#       # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
#       # 6th col variance of DL   "sig2_DL" 
#       # 7th col likelihood
#       colnames(NewY_mix_list)  <-c("curDur", "W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
#       NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_P_G)
#       NewY_mix_list$W_Ds = (1-NewY_mix_list$W_DL_PL) * (1-NewY_mix_list$W_P_G)
#       
#       # 1st col weight of global prior in De  "W_P_G"
#       # 2nd col mu_r # 3th col sig_r,  #4th col predY
#       # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
#       # 6th col variance of DL   "sig2_DL" 
#       # 7th col likelihood
#       colnames(pred_y_mix)  <-c("W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
#       pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_P_G)  #weight local prior in De
#       pred_y_mix$W_Ds = (1-pred_y_mix$W_DL_PL) * (1-pred_y_mix$W_P_G)  #weight of Ds in De
#       
#     }else if(modelname == 'PIM'){
#       NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
#       colnames(NewY_mix_list)  <-c("curDur","W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
#       # weight of global prior in De
#       NewY_mix_list$W_P_G = (1-NewY_mix_list$W_Ds) * (1-NewY_mix_list$W_PI_PL)
#       # weight of local prior in De
#       NewY_mix_list$W_L = (1- NewY_mix_list$W_Ds) * NewY_mix_list$W_PI_PL
#       
#       
#       # 1st col weight of Ds in De  "W_Ds"
#       # 2nd col mu_r # 3th col sig_r,  #4th col predY
#       # 5th col weight of local prior PL in integrated prior PI  "W_PI_PL" 
#       # 6th col variance of PI   "sig2_PI" 
#       # 7th col likelihood
#       colnames(pred_y_mix)  <-c("W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
#       # weight of global prior in De
#       pred_y_mix$W_P_G = (1-pred_y_mix$W_Ds) * (1-pred_y_mix$W_PI_PL)
#       # weight of local prior in De
#       pred_y_mix$W_L = (1- pred_y_mix$W_Ds) * pred_y_mix$W_PI_PL
#       
#     }
#   }
#   
#   PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model')]
#   PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
#   
#   
#   
#   Baypar = data.frame(
#     Nsub = subNo,
#     Exp = expName,
#     model = modelname,
#     phi = phi,
#     mu_p_s_log = mu_p_s_log,
#     mu_p_l_log = mu_p_l_log, 
#     mu_p_g = mu_p_g,
#     sig_pr2_s_log = sig_pr2_s_log,
#     sig_pr2_l_log = sig_pr2_l_log,
#     sig2_p_g = sig2_p_g,
#     sig2_mn = BR_par$sig2_mn,
#     sig2_t = BR_par$sig2_t,
#     looic = loo_1$looic,
#     p_loo = loo_1$p_loo,
#     elpd_loo = loo_1$elpd_loo,
#     se_looic = loo_1$se_looic,
#     se_p_loo = loo_1$se_p_loo,
#     waic = waic$waic,
#     p_waic =waic$p_waic,
#     se_waic = waic$se_waic,
#     se_p_waic = waic$se_p_waic,
#     elpd_waic = waic$elpd_waic
#   )
#   
#   Baypar_BR = data.frame(
#     Nsub = subNo,
#     Exp = expName,
#     model = modelname,
#     mu_p_s_log = mu_p_s_log,
#     mu_p_l_log = mu_p_l_log, 
#     sig_pr2_s_log = sig_pr2_s_log,
#     sig_pr2_l_log = sig_pr2_l_log,
#     sig2_mn = BR_par$sig2_mn,
#     sig2_t = BR_par$sig2_t,
#     mu_p_s_log_BR =  mu_p_s_log_BR,
#     mu_p_l_log_BR =  mu_p_l_log_BR,
#     sig_pr2_s_log_BR = sig_pr2_s_log_BR,
#     sig_pr2_l_log_BR = sig_pr2_l_log_BR,
#     mu_p_s_log_IR =  mu_p_s_log_IR,
#     mu_p_l_log_IR =  mu_p_l_log_IR,
#     sig_pr2_s_log_IR = sig_pr2_s_log_IR,
#     sig_pr2_l_log_IR = sig_pr2_l_log_IR
#   )
#   
#   if(modelname %in% c('DIM', 'LGM', 'PIM', 'IP')){
#     PredY_mix_list$W_P_S = PredY_mix_list$W_L   #weight local prior in De 
#     PredY_mix_list$W_P_L = PredY_mix_list$W_L
#     PredY_mix_list[which(PredY_mix_list$curDur <= 1),"W_P_L"] = 0
#     PredY_mix_list[which(PredY_mix_list$curDur >= 1),"W_P_S"] = 0
#     
#     NewY_mix_list$W_P_S = NewY_mix_list$W_L
#     NewY_mix_list$W_P_L = NewY_mix_list$W_L
#     NewY_mix_list[which(NewY_mix_list$curDur <= 1),"W_P_L"] = 0
#     NewY_mix_list[which(NewY_mix_list$curDur >= 1),"W_P_S"] = 0
#   }
#   NewY_mix_list$NSub =subNo
#   NewY_mix_list$exp = expName
#   NewY_mix_list$model = modelname
#   NewY_mix_list$group = 'mixed'
#   NewY_mix_list$part= partNo
#   PredY_mix_list$part= partNo
#   
#   return(list("Baypar" = Baypar, "NewY_mix_list" = NewY_mix_list, "PredY_mix_list" = PredY_mix_list, 
#               "mp" = mp, "phi_list" = phi_list, "Baypar_BR" = Baypar_BR))
# }



## The definition of the function to run Rstan model with sub data
funFitStanSub <- function(data_mix, myrstanModel, partNo, linear){
  
  # ## The definition of the function to run Rstan model with sub data
  # funFitStanSub <- function(data_mix){
  #   model = "LGM"
  #   myrstanModel <- rstan::stan_model(file=paste0("~/rp_globalprior_lrz/RStanCode/", model, "_linear.stan")) 
  #   partNo= 'part1'
  #   linear = 'linear'
  #   
  #   # library(rstanarm)
  #   # library(rstan)
  #   # library(tidyverse)
  #   # library(dplyr)
  #   # library(loo)
  modelname <-  unique(data_mix$model)
  subNo <- unique(data_mix$NSub)
  expName <- unique(data_mix$Exp)
  
  AllDat_Bayparlist_BR<- read.csv("RSTANMODELS/models/Baseline_log/AllDat_Bayparlist_BR.csv")
  if(linear == 'linear'){
    AllDat_Bayparlist_BR<- read.csv("RSTANMODELS/models/Baseline_linear/AllDat_Bayparlist_BR.csv")
  }
  BR_par = AllDat_Bayparlist_BR%>% filter(Nsub == subNo, Exp ==expName)
  print(paste0('Start run rstan model ', modelname,' on mixed session of Subject No.',subNo, ' in ', expName, ':', partNo ))
  
  durList <- sort(unique(data_mix$targetDur))
  for(i in 1: length(durList)){
    data_mix[which(data_mix$targetDur == durList[i]),"DurIdx"] = i
  }
  
  n_mix = length(data_mix$RP)
  xsnew <- seq(0.4, 0.8, 0.01)  #41
  xlnew <- seq(1.2, 2.4,0.01)  #121
  xnew <- c(seq(0.4, 0.8,0.01), seq(1.2, 2.4,0.01))   #162
  
  stan_data2 = NULL

  if(linear == 'linear'){
    stan_data2 = list( mu_p_s = BR_par$mu_p_s,
                       mu_p_l = BR_par$mu_p_l,
                       sig_pr2_s = BR_par$sig_pr2_s,
                       sig_pr2_l = BR_par$sig_pr2_l,
                       sig2_mn = BR_par$sig2_mn,
                       sig_t = BR_par$sig_t,
                       Y_mix = data_mix$RP,
                       n_mix = n_mix,
                       X_mix = data_mix$targetDur,
                       X_idx = data_mix$DurIdx,
                       xnew = xnew
    )  #data passed to stan
  }else{
    stan_data2 = list( mu_p_s_log = BR_par$mu_p_s_log,
                       mu_p_l_log = BR_par$mu_p_l_log, 
                       sig_pr2_s_log = BR_par$sig_pr2_s_log,
                       sig_pr2_l_log = BR_par$sig_pr2_l_log,
                       sig2_mn = BR_par$sig2_mn,
                       sig2_t = BR_par$sig2_t,
                       Y_mix = data_mix$RP, 
                       n_mix = n_mix, 
                       X_mix = data_mix$targetDur,
                       X_idx = data_mix$DurIdx,
                       xnew = xnew)  #data passed to stan
  }
  
  parameters2 <- c("mu_p_g", "sig2_p_g", "ynew_mix", "ypred_mix")
  myinits2 <- list(
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1),    
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1),
    list(mu_p_g=1, sig2_p_g = 0.1)
  )
  
  
  if (modelname == 'BPM'){
    parameters2 <- c("mp", "ynew_mix") 
    myinits2 <- list(list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5),
                     list(mp = 0.5))
  }
  
  if (modelname  == 'DIM'){
    parameters2 <- c( "mu_p_g", "sig2_p_g", "phi", "ynew_mix", "ypred_mix") 
  }
  
  if (modelname == 'IP') {
    parameters2 <- c("mu_p_s", "mu_p_l","sig_pr2_s", "sig_pr2_l", "ynew_mix", "ypred_mix") 
    myinits2 <- list(
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),    
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1),
      list(mu_p_s=0.7, sig_pr2_s = 0.1, mu_p_l=1.2, sig_pr2_l =0.1)
    )
    
    stan_data2 = list(sig2_mn = BR_par$sig2_mn,
                      sig2_t = BR_par$sig2_t,
                      Y_mix = data_mix$RP,
                      n_mix = n_mix,
                      X_mix = data_mix$targetDur,
                      xnew = xnew
    )  #data passed to stan
    
    if(linear == 'linear'){
      stan_data2 = list(sig2_mn = BR_par$sig2_mn,
                        sig_t = BR_par$sig_t,
                        Y_mix = data_mix$RP,
                        n_mix = n_mix,
                        X_mix = data_mix$targetDur,
                        xnew = xnew
      )  #data passed to stan
    }
    
  }
  
  
  subfit2 <- sampling(myrstanModel, 
                      stan_data2,
                      init=myinits2,
                      iter=12000,
                      chains=8,
                      thin=1,
                      control = list(adapt_delta = 0.9,
                                     max_treedepth = 15),
                      cores = parallel::detectCores())
  
  fitpar <- summary(subfit2, pars = parameters2)$summary
  list_of_draws <- rstan::extract(subfit2, pars = parameters2)
  phi = 0
  mp = NULL
  phi_list = NULL
  Baypar = NULL
  NewY_mix_list = NULL
  PredY_mix_list = NULL
  
  log_lik_rlt <- extract_log_lik(subfit2)
  loo_1 <- loo(log_lik_rlt)
  waic = waic(log_lik_rlt)
  mu_p_g =  0
  sig2_p_g = 0  
  mu_p_s_IR =  0
  mu_p_l_IR =  0
  sig_pr2_s_IR = 0
  sig_pr2_l_IR = 0
  if(modelname == 'IP'){
    mu_p_s_IR =  mean(list_of_draws$mu_p_s)
    mu_p_l_IR =  mean(list_of_draws$mu_p_l)
    sig_pr2_s_IR = mean(list_of_draws$sig_pr2_s)
    sig_pr2_l_IR = mean(list_of_draws$sig_pr2_l)
  }else{
    if(modelname != 'BPM'){
      mu_p_g =  mean(list_of_draws$mu_p_g)
      sig2_p_g = mean(list_of_draws$sig2_p_g)
    }
  }
  
  
  if(modelname == 'BPM'){
    mp0 =  mean(list_of_draws$mp)
    mp <-  c(mp0, 1-mp0, 0, subNo, expName, modelname)
    mp = data.frame(mp)
    
    pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
    y_pred_mix_list <- list_of_draws$ypred_mix
    for (i in 1:7){
      for (j in 1:n_mix){
        pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
      }
    }
    pred_y_mix = data.frame(pred_y_mix)
    colnames(pred_y_mix)  <-c("predY_s", "predY_l", "predY", "mu_r", "wp_ps", "wp_pl", "log_lik")
    
    PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model', 'DurIdx')]
    PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
    
    PredY_mix_list$mp = mp0
    PredY_mix_list$W_P_S = mp0*PredY_mix_list$wp_ps
    PredY_mix_list$W_P_L = (1-mp0)* PredY_mix_list$wp_pl
    PredY_mix_list$W_Ds = mp0*(1-PredY_mix_list$wp_ps)+ (1-mp0)*(1-PredY_mix_list$wp_pl)
    PredY_mix_list$W_P_G = 0
    PredY_mix_list$group ='mixed'
    NewY_mix_list$part= partNo
    PredY_mix_list$sig_r =0. # Question: how to calculate sig_r in BPM model?
  }
  
  if(modelname == 'DIM'){
    phi = mean(list_of_draws$phi)
    ynew_mix_list <- list_of_draws$ynew_mix
    newy_mix <- matrix(rep(0, 162, 9), nrow = 162, ncol = 9)
    for (i in 1:9){
      for (j in 1:162){
        newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
      }
    }
    
    NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
    # 1st col //weight of D_G in De  "W_DG"
    # 2nd col mu_r # 3th col sig_r,  #4th col predY
    # 5th col weight of local prior PL in DL  "W_DL_PL"  
    # 6th col weight of global prior in DG  "W_DG_PG"
    # 7th col variance of DL   "sig2_DL" 
    # 8th col variance of DG  "sig2_DG"
    # 9th col likelihood
    colnames(NewY_mix_list)  <-c("curDur","W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
    #weight of global prior in De
    NewY_mix_list$W_P_G =  NewY_mix_list$W_DG_PG * NewY_mix_list$W_DG
    #weight of local prior in De
    NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_DG)
    #weight of Ds in De 
    NewY_mix_list$W_Ds =(1-NewY_mix_list$W_DL_PL) *(1-NewY_mix_list$W_DG) + NewY_mix_list$W_DG* (1- NewY_mix_list$W_DG_PG)
    
    pred_y_mix <- matrix(rep(0, n_mix, 9), nrow = n_mix, ncol = 9)
    y_pred_mix_list <- list_of_draws$ypred_mix
    for (i in 1:9){
      for (j in 1:n_mix){
        pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
      }
    }
    pred_y_mix = data.frame(pred_y_mix)
    # 1st col weight of D_G in De  "W_DG"
    # 2nd col mu_r # 3th col sig_r,  #4th col predY
    # 5th col weight of local prior PL in DL  "W_DL_PL"  
    # 6th col weight of global prior in DG  "W_DG_PG"
    # 7th col variance of DL   "sig2_DL" 
    # 8th col variance of DG  "sig2_DG"
    # 9th col likelihood
    colnames(pred_y_mix)  <-c("W_DG", "mu_r", "sig_r", "predY", "W_DL_PL","W_DG_PG", "sig2_DL","sig2_DG","log_lik")
    #weight of global prior in De
    pred_y_mix$W_P_G =  pred_y_mix$W_DG_PG * pred_y_mix$W_DG
    #weight of local prior in De
    pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_DG)
    #weight of Ds in De 
    pred_y_mix$W_Ds = (1-pred_y_mix$W_DL_PL) *(1-pred_y_mix$W_DG) + pred_y_mix$W_DG* (1-pred_y_mix$W_DG_PG)
    
  }
  
  if(modelname %in% c('LGM', 'PIM', 'IP')){
    ynew_mix_list <- list_of_draws$ynew_mix
    newy_mix <- matrix(rep(0, 162, 7), nrow = 162, ncol = 7)
    for (i in 1:7){
      for (j in 1:162){
        newy_mix[j,i] <-  mean(ynew_mix_list[,j,i] )
      }
    }
    
    pred_y_mix <- matrix(rep(0, n_mix, 7), nrow = n_mix, ncol = 7)
    y_pred_mix_list <- list_of_draws$ypred_mix
    for (i in 1:7){
      for (j in 1:n_mix){
        pred_y_mix[j,i] <-  mean(y_pred_mix_list[,j,i] )
      }
    }
    pred_y_mix = data.frame(pred_y_mix)
    
    if(modelname%in% c('LGM','IP')){
      NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
      # 1st col weight of global prior in De  "W_P_G"
      # 2nd col mu_r # 3th col sig_r,  #4th col predY
      # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
      # 6th col variance of DL   "sig2_DL" 
      # 7th col likelihood
      colnames(NewY_mix_list)  <-c("curDur", "W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
      NewY_mix_list$W_L = NewY_mix_list$W_DL_PL * (1-NewY_mix_list$W_P_G)
      NewY_mix_list$W_Ds = (1-NewY_mix_list$W_DL_PL) * (1-NewY_mix_list$W_P_G)
      
      # 1st col weight of global prior in De  "W_P_G"
      # 2nd col mu_r # 3th col sig_r,  #4th col predY
      # 5th col weight of local prior PL in DL  "W_DL_PL" = 1-"W_DL_Ds"
      # 6th col variance of DL   "sig2_DL" 
      # 7th col likelihood
      colnames(pred_y_mix)  <-c("W_P_G", "mu_r", "sig_r", "predY", "W_DL_PL","sig2_DL","log_lik")
      pred_y_mix$W_L = pred_y_mix$W_DL_PL * (1-pred_y_mix$W_P_G)  #weight local prior in De
      pred_y_mix$W_Ds = (1-pred_y_mix$W_DL_PL) * (1-pred_y_mix$W_P_G)  #weight of Ds in De
      
    }else if(modelname == 'PIM'){
      NewY_mix_list <- data.frame(cbind(xnew, newy_mix))
      colnames(NewY_mix_list)  <-c("curDur","W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
      # weight of global prior in De
      NewY_mix_list$W_P_G = (1-NewY_mix_list$W_Ds) * (1-NewY_mix_list$W_PI_PL)
      # weight of local prior in De
      NewY_mix_list$W_L = (1- NewY_mix_list$W_Ds) * NewY_mix_list$W_PI_PL

      # 1st col weight of Ds in De  "W_Ds"
      # 2nd col mu_r # 3th col sig_r,  #4th col predY
      # 5th col weight of local prior PL in integrated prior PI  "W_PI_PL" 
      # 6th col variance of PI   "sig2_PI" 
      # 7th col likelihood
      colnames(pred_y_mix)  <-c("W_Ds", "mu_r", "sig_r", "predY", "W_PI_PL","sig2_PI", "log_lik")
      # weight of global prior in De
      pred_y_mix$W_P_G = (1-pred_y_mix$W_Ds) * (1-pred_y_mix$W_PI_PL)
      # weight of local prior in De
      pred_y_mix$W_L = (1- pred_y_mix$W_Ds) * pred_y_mix$W_PI_PL
    }
  }
  
  PredY_mix_list <- data_mix[c('id','trialnum','NSub','targetDur', 'RP','Exp','group', 'model')]
  PredY_mix_list <- cbind(PredY_mix_list, pred_y_mix)
  
  
  
  Baypar = data.frame(
    Nsub = subNo,
    Exp = expName,
    model = modelname,
    ver = linear,
    phi = phi,
    mu_p_g = mu_p_g,
    sig2_p_g = sig2_p_g,
    looic = loo_1$looic,
    p_loo = loo_1$p_loo,
    elpd_loo = loo_1$elpd_loo,
    se_looic = loo_1$se_looic,
    se_p_loo = loo_1$se_p_loo,
    waic = waic$waic,
    p_waic =waic$p_waic,
    se_waic = waic$se_waic,
    se_p_waic = waic$se_p_waic,
    elpd_waic = waic$elpd_waic,
    mu_p_s_IR =  mu_p_s_IR,
    mu_p_l_IR =  mu_p_l_IR,
    sig_pr2_s_IR = sig_pr2_s_IR,
    sig_pr2_l_IR = sig_pr2_l_IR,
    part = partNo
  )
  
  
  if(modelname %in% c('DIM', 'LGM', 'PIM', 'IP')){
    PredY_mix_list$W_P_S = PredY_mix_list$W_L   #weight local prior in De 
    PredY_mix_list$W_P_L = PredY_mix_list$W_L
    PredY_mix_list[which(PredY_mix_list$curDur <= 1),"W_P_L"] = 0
    PredY_mix_list[which(PredY_mix_list$curDur >= 1),"W_P_S"] = 0
    
    NewY_mix_list$W_P_S = NewY_mix_list$W_L
    NewY_mix_list$W_P_L = NewY_mix_list$W_L
    NewY_mix_list[which(NewY_mix_list$curDur <= 1),"W_P_L"] = 0
    NewY_mix_list[which(NewY_mix_list$curDur >= 1),"W_P_S"] = 0
  }
  NewY_mix_list$NSub =subNo
  NewY_mix_list$exp = expName
  NewY_mix_list$model = modelname
  NewY_mix_list$group = 'mixed'
  NewY_mix_list$part= partNo
  PredY_mix_list$part= partNo
  
  return(list("Baypar" = Baypar, "NewY_mix_list" = NewY_mix_list, "PredY_mix_list" = PredY_mix_list, 
              "mp" = mp, "phi_list" = phi_list))
}
