#detach(package:dplyr)
#detach(package:plyr)
#detach(package:neuralnet)
library(tidyverse)
library(ez)
library(nleqslv)
library(ggsignif) 
library(ggpubr)
library(R.matlab)
library(rstatix)
library(reshape2)
library(latex2exp)
library(dplyr)
library(rlist)
library(loo)
usePackage <- function(pk){
  for (p in pk) {
    if (!is.element(p, installed.packages()[,1])) install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}
usePackage(c('knitr','R.matlab','data.table','xtable','dplyr','tidyr',
             'grid','gridExtra','broom','ggplot2','nlme','ez','RColorBrewer',
             'scales','ggsignif','cowplot','tidyverse','ggsignif','apaTables'))


theme_set(theme_classic())

# set global options
mytheme <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"))+
  theme_classic() + 
  theme(strip.background = element_blank()) 


## calulate geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#------------------merger and load valid data--------------
# define the function to read all csv files after Preprocessing and merge them together 
myMergeData <- function(cvsfiles, filename){
  merge.data = read.csv(file.path(cvsfiles[1]), header=T,sep=",")
  if (length(cvsfiles) >= 2) {
    for (i in 2:length(cvsfiles)){
      new.data = read.csv(file.path(cvsfiles[i]), header=T, sep=",")
      merge.data = rbind(merge.data,new.data)
    }
  }
  filename1 <- paste0("data/AllData_", filename, ".csv")
  write.csv(file=filename1, merge.data)
}

RawdataPath <- "data/RawData/";

#Merge all of the Trials data
AllDataMerge <- function(files){
  merge.data = read.csv(file.path(paste0("data/AllData_", files[1], ".csv")), header=T, sep=",")
  if (length(files) >= 2) {
    for (i in 2:length(files)){
      new.data = read.csv(file.path(paste0("data/AllData_", files[i], ".csv")), header=T, sep=",")
      merge.data = rbind(merge.data,new.data)
    }
  }
  filename <- paste0("data/AllExpData.csv")
  write.csv(file=filename, merge.data)
}



Correct_Colnames <- function(df) {
  delete.columns <- grep("X|(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
  if (length(delete.columns) > 0) {
    #row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
    df <- df[,-delete.columns]
    colnames(df) <- gsub("^X", "",  colnames(df))
  }
  return(df)
}

## read in data, factorize group subject
readAllData <- function(filepath, rm_outlier, priortype){
  expdata = read.csv(filepath)
  expdata <- Correct_Colnames(expdata)
  
  expdata <- expdata %>% 
    dplyr::filter(abs(expdata$phyTargetDur - expdata$targetDur) <= 30)%>% 
    dplyr::filter(expdata$valid==1) 
  
  expdata$firstSession = factor(expdata$order, levels = c(1, 2), labels = c('short','long'))
  
  #add two new columns 'priortype' and 'range' isntead of 'group'
  expdata$group = factor(expdata$group, labels = c('short','long', 'mixed'))
  expdata$priortype = 'BR'
  expdata[which(expdata$group == 'mixed'),"priortype"] = 'IR'
  expdata$priortype = as.factor(expdata$priortype)
  expdata$range = TRUE
  expdata[which(expdata$targetDur > 1),"range"] = FALSE
  expdata$range = as.factor(expdata$range)
  
  dat_exps = expdata
  if(rm_outlier){
    dat_exps = NULL
    exp_outliers<- expdata %>% dplyr::group_by(Exp, targetDur, range, priortype)%>%
      dplyr::summarise(limit_low = quantile(repError, 0.025),  #0.025
                       limit_high = quantile(repError, 0.975))  #0.975
    expdata <- left_join(expdata, exp_outliers, by = c("Exp", "targetDur", "range", "priortype") )
    exp_outliers_sub<- expdata %>% dplyr::group_by(NSub, Exp, targetDur, range, priortype)%>%
      dplyr::summarise(limit_low_sub = quantile(repError, 0.025),  #0.025
                       limit_high_sub = quantile(repError, 0.975))  #0.975
    expdata <- left_join(expdata, exp_outliers_sub, by = c("NSub", "Exp", "targetDur", "range", "priortype") )
   
     dat_exps <- expdata %>% filter( repError < limit_high, 
                                     repError > limit_low,
                                     repError < limit_high_sub, 
                                     repError > limit_low_sub)
     
     ratio <- length(dat_exps$NB)/length(expdata$NB) #0.8957678approximately  10% of all trials
  }

if(priortype){
  dat_exps$range = 'short'
  dat_exps[which(dat_exps$targetDur > 1),"range"] = 'long'
  dat_exps$range = as.factor(dat_exps$range)
}else{
  dat_exps$range = NULL
  dat_exps$priortype = NULL
}
return(dat_exps)
}



#convert matlab data to csv files
convertMatToCSV<- function(matfiles, foldername){
  csvfiles <- {}
  for (i in 1:length(matfiles)){
    rawDat <- readMat(file.path(getwd(),  paste0(RawdataPath,foldername,"/", matfiles[i])))
    newfilename <- paste0(substring(matfiles[i],1,which(strsplit(matfiles[i], "")[[1]]==".")), 'csv')
    csvfilepath <- paste0("data/trials/",foldername,"/", newfilename)
    write.csv(rawDat$trials, file = csvfilepath)
    csvfiles= rbind2(csvfiles, newfilename);
  }
  return(csvfiles);
}


# summarize the subject infomation for each expeirment
summariseSubInfo <- function(matfiles, foldername){
  expInfoList <- data.frame()
  csvfilepath <- paste0("data/subInfo_", foldername, '.csv')
  for (i in 1:length(matfiles)){
    rawDat <- readMat(file.path(getwd(),  paste0(RawdataPath,foldername,"/", matfiles[i])))
    expInfo <- rawDat$expInfo
    expPar <- rawDat$aData
    expInfoList <- rbind2(expInfoList, c(expInfo[2], expInfo[3], expPar[1], expPar[2]));
  }
  colnames(expInfoList) <- c("Age", "Gender", "Order", "Loc") #'First: S(1) or L(2)', 'Loc: S(1) or L(2)'};
  expInfoList$Age = as.numeric(expInfoList$Age)
  expInfoList$Name = paste0('Sub', 1:length(matfiles))
  expInfoList <-data.frame(expInfoList)
  write.csv(expInfoList, file = csvfilepath)
  return(expInfoList);
}

#summerize the subject information for all experient together
summariseSubInfoAll <-function(ExpNames){
  merge.data = read.csv(file.path(paste0("data/subInfo_", ExpNames[1], ".csv")), header=T, sep=",")
  merge.data$Exp <- ExpNames[1]
  if (length(ExpNames) >= 2) {
    for (i in 2:length(ExpNames)){
      new.data = read.csv(file.path(paste0("data/subInfo_", ExpNames[i], ".csv")), header=T, sep=",")
      new.data$Exp <- ExpNames[i]
      merge.data = rbind(merge.data,new.data)
    }
  }
  filename <- paste0("data/AllExpSubInfo.csv")
  write.csv(file=filename, merge.data)
}


#summarize the data for each subject group by conditions gp_var
summarizedata<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  sdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarise(m_err = mean(repError),
                     m_RP = mean(RP), 
                     m_BIAS = mean(BIAS),  # abs(repError)
                     n = n(), 
                     se_err = sd(repError)/ sqrt(n-1),
                     sd_err = sd(repError),
                     se_RP = sd(RP)/ sqrt(n-1),
                     sd_RP = sd(RP),
                     CV = sd_RP/m_RP, 
                     NBIAS = m_BIAS/mean(targetDur),
                     NRMSE = sqrt(NBIAS^2+CV^2))
  return(sdata)
}


#summarize the predicated RP group by conditions gp_var
summarizePredY<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  sdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarize(n = n(),
                     m_m_RP_BIAS = mean(m_RP_BIAS),
                     m_m_prederr = mean(m_prederr),
                     m_m_predBIAS = mean(m_predBIAS),  #prediction minus target duration
                     m_m_predNBIAS =mean(m_predNBIAS),
                     m_m_RP = mean(m_RP),
                     m_cv_observed = mean(cv_observed),
                     m_sd_RP = mean(sd_RP),
                     mm_mu_r = mean(m_mu_r),
                     mm_sig_r = mean(m_sig_r),
                     mm_pred_Var = mean(m_sig_r- sd_RP),
                     m_cv_Pred = mean(cv_Pred),
                     se_prederr = sd(m_prederr)/sqrt(n-1),
                     se_pred_Var = sd(m_sig_r- sd_RP)/sqrt(n-1),
                     predRP_err = mean(m_mu_r-m_RP),
                     predVar_err = mean(m_sig_r-sd_RP),
                     predRP_rerr = mean(abs(m_mu_r-m_RP)/m_RP),
                     predVar_rerr = mean(abs(m_sig_r-sd_RP)/sd_RP),
                     pred_cv = mean(m_sig_r/m_mu_r),
                     predcv_err = pred_cv-m_cv_observed,
                     predcv_rerr = mean(abs(pred_cv-m_cv_observed)/m_cv_observed)
                     )
  return(sdata)
}

#summarize the predicated RP group by conditions gp_var
summarizePredY2<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  sdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarize(n = n(),
                     m_RP = mean(RP), 
                     sd_RP = sd(RP),
                     m_mu_r = mean(mu_r),
                     m_sig_r = mean(sig_r),
                     log_lik =mean(log_lik),
                     cv = sd_RP/m_RP,
                     pred_cv = mean(sig_r/mu_r),
                     predRP_err = mean(mu_r-RP),
                     se_predRP_err = sd(mu_r-RP)/sqrt(n-1),
                     predVar_err = mean(m_sig_r-sd_RP),
                     predRP_rerr = mean(abs(m_mu_r-m_RP)/m_RP),
                     predVar_rerr = mean(abs(m_sig_r-sd_RP)/sd_RP),
                     predcv_err = pred_cv-cv,
                     predcv_rerr = mean(abs(pred_cv-cv)/cv))
  
  return(sdata)
}

summarizeNewY2<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  sdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarize(n = n(),
                     m_mu_r = mean(mu_r),
                     m_sig_r = mean(sig_r),
                     log_lik =mean(log_lik))
  
  return(sdata)
}



#summarize the data for each subject group by conditions gp_var
summarizeMeanData<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  sdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarise(m_err = mean(repError),
                     m_RP = mean(RP), 
                     m_BIAS = mean(BIAS),
                     n = n(), 
                     se_err = sd(repError)/ sqrt(n-1),
                     sd_err = sd(repError),
                     se_BIAS = sd(BIAS)/sqrt(n-1),
                     sd_BIAS = sd(repError),
                     se_RP = sd(RP)/ sqrt(n-1),
                     sd_RP = sd(RP),
                     CV = sd_RP/m_RP, 
                     NBIAS = m_BIAS/mean(targetDur),
                     NRMSE = sqrt(NBIAS^2+CV^2))%>%
    dplyr::group_by(!!!gp_vars[2:length(gp_vars)]) %>%
    dplyr::summarise( m_m_BIAS = mean(m_BIAS),
                      m_m_err = mean(m_err),
                      m_m_RP = mean(m_RP),
                      m_sd_RP = mean(sd_RP),
                      m_sd_err = mean(sd_err),
                      m_se_RP = mean(se_RP),
                      m_se_err = mean(se_err),
                      m_sd_BIAS = mean(sd_BIAS),
                      m_se_BIAS = mean(se_BIAS),
                      m_CV = mean(CV),
                      n= n(),
                      se_CV = sd(CV)/sqrt(n-1),
                      m_NBIAS = mean(NBIAS),
                      m_NRMSE = mean(NRMSE))
  return(sdata)
}



#summarize the data for each subject group by conditions gp_var
mm_summarizedata<- function(df, gp_var) {
  gp_vars = syms(gp_var)
  mdata <- df %>%
    dplyr::group_by(!!!gp_vars) %>%
    dplyr::summarise(m_err = mean(repError),
                     m_RP = mean(RP), 
                     m_BIAS = mean(BIAS),
                     n = n(), 
                     se_err = sd(repError)/ sqrt(n-1),
                     sd_err = sd(repError),
                     se_BIAS = sd(BIAS)/sqrt(n-1),
                     sd_BIAS = sd(repError),
                     se_RP = sd(RP)/ sqrt(n-1),
                     sd_RP = sd(RP),
                     CV = sd_RP/m_RP, 
                     NBIAS = m_BIAS/mean(targetDur),
                     NRMSE = sqrt(NBIAS^2+CV^2))%>%
    dplyr::group_by(!!!gp_vars[2:length(gp_vars)]) %>%
    dplyr::summarise( m_m_BIAS = mean(m_BIAS),
                      m_m_err = mean(m_err),
                      m_m_RP = mean(m_RP),
                      m_sd_RP = mean(sd_RP),
                      m_sd_err = mean(sd_err),
                      m_se_RP = mean(se_RP),
                      m_se_err = mean(se_err),
                      m_sd_BIAS = mean(sd_BIAS),
                      m_se_BIAS = mean(se_BIAS),
                      m_CV = mean(CV),
                      n= n(),
                      se_CV = sd(CV)/sqrt(n-1),
                      m_NBIAS = mean(NBIAS),
                      m_NRMSE = mean(NRMSE))%>%
    dplyr::group_by(!!!gp_vars[3:length(gp_var)]) %>%
    dplyr::summarise(m_m_BIAS = mean(m_m_BIAS),
                     m_m_err = mean(abs(m_m_err)),
                     m_m_RP = mean(m_m_RP),
                     m_sd_RP = mean(m_sd_RP),
                     m_sd_err = mean(m_sd_err),
                     m_se_RP = mean(m_se_RP),
                     m_se_BIAS = mean(m_se_BIAS),
                     m_NBIAS = mean(m_NBIAS),
                     m_CV = mean(m_CV),
                     se_CV = mean(se_CV),
                     m_NRMSE = mean(m_NRMSE))
  return(mdata)
}



mRep_model <- function(df) {
  lm(RP ~ targetDur, data = df)
}


mRP_model <- function(df) {
  lm(m_m_RP ~ targetDur, data = df)
}

# definition of fuction calculate slope
getSlope <- function(df, gp_var){
  gp_vars = syms(gp_var) 
  slopes <- df %>% 
    dplyr::group_by(!!!gp_vars) %>% nest()  %>%  # nested data
    mutate(model = map(data, mRep_model)) %>%  # linear regression
    mutate(slope = map(model, broom::tidy)) %>%  # get estimates
    unnest(slope, .drop = TRUE) %>% # remove raw data
    select(-std.error,-statistic, -p.value) %>%  # remove unnessary columns
    spread(term, estimate) %>%   # spread stimates
    dplyr::rename(minRP = `(Intercept)`, slope = targetDur)  # rename columns
  slopes$data = NULL
  slopes$model = NULL
  slopes$inP = slopes$minRP/(1-slopes$slope)
  
  # get the max and min of target duration
  dat_maxDur =  df %>% 
    dplyr::group_by(Exp) %>%
    dplyr::summarise(maxDur = max(targetDur),  mixDur = min(targetDur))
  
  slopes <-  left_join(slopes, dat_maxDur, by = c('Exp'))
  
  # find out the indifferent point outliers
  slopes$inPOutlier = FALSE
  slopes[which(slopes$inP > slopes$maxDur |slopes$inP < slopes$mixDur),"inPOutlier"] = TRUE
  slopes$inPOutlier =as.factor(slopes$inPOutlier)
  return(slopes)
}

# definition of fuction to calculate slope
getmSlope <- function(df, gp_var){
  mslope <- getSlope(df, gp_var)
  mslope$CI <- 1- mslope$slope 
  slope <- mslope %>%# filter(inPOutlier == FALSE)%>%  #exclude the outliers
    group_by(!!!syms(gp_var[2:length(gp_var)])) %>% 
    dplyr::summarise(m_slope = mean(slope),
                     m_CI = mean(CI),
                     m_minRP = mean(minRP),
                     m_inP = mean(inP),
                     n = n(), 
                     se_slope = sd(slope)/ sqrt(n-1),
                     se_CI = sd(CI)/sqrt(n-1)
    )
  return(slope)
}


getmmSlope <- function(df, gp_var){
  gp_vars = syms(gp_var) 
  mm_slope <- df %>% 
    dplyr::group_by(!!!gp_vars) %>% nest()  %>%  # nested data
    mutate(model = map(data, mRP_model)) %>%  # linear regression
    mutate(slope = map(model, broom::tidy)) %>%  # get estimates out
    unnest(slope, .drop = TRUE) %>% # remove raw data
    select(-std.error,-statistic, -p.value) %>%  # remove unnessary clumns
    spread(term, estimate) %>%   # spread stimates
    dplyr::rename(minRP = `(Intercept)`, slope = targetDur)  # rename columns
  mm_slope$inP = mm_slope$minRP /(1-mm_slope$slope)
  
  # get the max and min of target duration
  dat_maxDur =  df %>% 
    dplyr::group_by(Exp) %>%
    dplyr::summarise(maxDur = max(targetDur),  mixDur = min(targetDur))
  
  mm_slope <-  left_join(mm_slope, dat_maxDur, by = c('Exp'))
  
  # find out the indifferent point outliers
  mm_slope$inPOutlier = FALSE
  mm_slope[which(mm_slope$inP > mm_slope$maxDur |mm_slope$inP < mm_slope$mixDur),"inPOutlier"] = TRUE
  mm_slope$inPOutlier =as.factor(mm_slope$inPOutlier)
  return(mm_slope)
}




# Custom function to find indifference point 
library('boot')
getInP <- function(df, idx){
  vars <- c('NSub', 'Exp', 'priortype', 'range')
  gp_vars = syms(vars) 
  slopes <- df[idx, ] %>% 
    dplyr::group_by(!!!gp_vars) %>% nest()  %>%  # nested data
    mutate(model = map(data, mRep_model)) %>%  # linear regression
    mutate(slope = map(model, broom::tidy)) %>%  # get estimates out
    unnest(slope, .drop = TRUE) %>% # remove raw data
    select(-std.error,-statistic, -p.value) %>%  # remove unnessary clumns
    spread(term, estimate) %>%   # spread stimates
    dplyr::rename(minRP = `(Intercept)`, slope = targetDur)  # rename columns
  slopes$inP = slopes$minRP /(1-slopes$slope)
  return(slopes$inP)
}




#Custom function to get slope 
getBootSlope <- function(df, idx){
  vars <- c('NSub', 'Exp', 'priortype', 'range')
  gp_vars = syms(vars) 
  slopes <- df[idx, ] %>% 
    dplyr::group_by(!!!gp_vars) %>% nest()  %>%  # nested data
    mutate(model = map(data, mRep_model)) %>%  # linear regression
    mutate(slope = map(model, broom::tidy)) %>%  # get estimates out
    unnest(slope, .drop = TRUE) %>% # remove raw data
    select(-std.error,-statistic, -p.value) %>%  # remove unnessary clumns
    spread(term, estimate) %>%   # spread stimates
    dplyr::rename(minRP = `(Intercept)`, slope = targetDur)  # rename columns
  return(slopes$slope)
}




#Custom function to get Intercept 
getIntercept <- function(df, idx){
  vars <- c('NSub', 'Exp', 'priortype', 'range')
  gp_vars = syms(vars) 
  slopes <- df[idx, ] %>% 
    dplyr::group_by(!!!gp_vars) %>% nest()  %>%  # nested data
    mutate(model = map(data, mRep_model)) %>%  # linear regression
    mutate(slope = map(model, broom::tidy)) %>%  # get estimates out
    unnest(slope, .drop = TRUE) %>% # remove raw data
    select(-std.error,-statistic, -p.value) %>%  # remove unnessary clumns
    spread(term, estimate) %>%   # spread stimates
    dplyr::rename(minRP = `(Intercept)`, slope = targetDur)  # rename columns
  return(slopes$minRP)
}
