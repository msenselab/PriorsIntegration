source('ana_functions.R')

#########Data Preprocessing#################

Exp2_mat = c(
  'AY3626M2.mat',
  'DZ4423F2.mat',
  'FJW4529M1.mat',
  'JGP3122M1.mat',
  'MK3426M1.mat',
  'MK3523M1.mat',
  'SA3027M2.mat',
  'TI3324F2.mat',
  'WM3221M2.mat',
  'WPQ23F2.mat',
  'YW4124M1.mat',
  'ZH4124F1.mat',
  'PI24M1.mat', 
  'MA26M2.mat',
  'DJ25M1.mat',
  'CV31F2.mat'
)


Exp1_mat = c(
  'BK26M2.mat',
  'CB25F1.mat',
  'DY4723F2.mat',
  'JLL4823F1.mat',
  'JP4930F1.mat',
  'PD23F1.mat',
  'RB36F1.mat',
  'SAS3734F1.mat',
  'SS27M2.mat',
  'XC21F2.mat',
  'XH24F1.mat',
  'YuchenChi22F1.mat',
  'ZY4622F2.mat', 
  'GV25F2.mat',
  'JK23M2.mat',
  'MK24F2.mat'
)




Exp3_mat = c(
  'SY0125F1.mat',
  'XZ0223F1.mat',
  'SA0327M2.mat',
  'JM0432M2.mat',
  'AH0526M1.mat',
  'PP0626M1.mat',
  'SC0719F2.mat',
  'SY0125F1.mat'
)


## defination of the function converting  trials to tables with additional columns
##parameters:
##  subfiles-- list of files
##  filepath -- subfoldername
##  blockTrials -- trials number of a block
##  pracBlocks -- the number of practise block

myDatPreprocessing <- function(subfiles, subInfo, filepath, ExpName, blockTrials, pracBlocks){
  ProjectDir <- getwd()
  mypath <- paste0("data/trials/",filepath)
  dataDir <- file.path(ProjectDir, mypath)
  cvsfiles ={}
  Location <- subInfo['Loc'][[1]]
  Order<- subInfo['Order'][[1]]  #1: short group as first session, 2: long group as first session
  for (i in 1:length(subfiles)){
 
    data_Trail <- data.frame(read.csv(file.path(dataDir, subfiles[i]), header=T, sep=",")) 
    dataTrial <- Correct_Colnames(data_Trail)
    colnames(data_Trail) <- c("trialnum","group", "interval", "targetDur", "phyTargetDur", "RP")
    data_Trail$NSub = i
    data_Trail$NT<-1:nrow(data_Trail)

    #set Exp
    data_Trail$Exp= ExpName
    
    #set NB
    data_Trail$NB= 0
    data_Trail[which(data_Trail$NT<= blockTrials*pracBlocks),"NB"] = 0
    
    #set repError 
    data_Trail$repError = data_Trail$RP - data_Trail$targetDur;
    
    #set valid 
    data_Trail$valid= 1
    data_Trail[which(data_Trail$NT< blockTrials * pracBlocks),"valid"] = 0
     dat <- data_Trail %>% dplyr::filter(data_Trail$valid == 1)
    data_Trail[which(data_Trail$NT<= blockTrials*pracBlocks),"NB"] = 0
    ntstart<- dat$NT[1]
    ntend <- dat$NT[length(dat$NT)]
    
    for (index in ntstart:ntend){
      data_Trail[which(data_Trail$NT==index),"NB"] = floor((index-ntstart) / blockTrials)+1 
    }
    
    data_Trail[which( abs(data_Trail$phyTargetDur- data_Trail$targetDur)> 20), "valid"] = 0
    data_Trail$firstLoc = as.factor(Location[i])
    data_Trail$order = as.factor(Order[i])
    filename = paste0("data/csv/", filepath, "/NSub", toString(i), "_", subfiles[i])
    cvsfiles= rbind2(cvsfiles, filename);
    write.csv(file=file.path(getwd(), filename), data_Trail)
  }
  return(cvsfiles);
}


#exp1
file_Exp1 <- 'Exp1'
subfiles_Exp1 <- convertMatToCSV(Exp1_mat, file_Exp1)
subInfo_Exp1 <- summariseSubInfo(Exp1_mat, file_Exp1)
cvsfiles_1 <- myDatPreprocessing(subfiles_Exp1, subInfo_Exp1, file_Exp1, "Exp1", 24, 1)


#exp2
file_Exp2 <- 'Exp2'
subfiles_Exp2 <- convertMatToCSV(Exp2_mat, file_Exp2)
subInfo_Exp2 <- summariseSubInfo(Exp2_mat, file_Exp2)
cvsfiles_2 <- myDatPreprocessing(subfiles_Exp2, subInfo_Exp2, file_Exp2, "Exp2", 24, 1)

#exp3
file_Exp3 <- 'Exp3'
subfiles_Exp3 <- convertMatToCSV(Exp3_mat, file_Exp3)
subInfo_Exp3 <- summariseSubInfo(Exp3_mat, file_Exp3)
cvsfiles_3 <- myDatPreprocessing(subfiles_Exp3, subInfo_Exp3, file_Exp3, "Exp3", 24, 1)


myMergeData(cvsfiles_1, file_Exp1)
myMergeData(cvsfiles_2, file_Exp2)
myMergeData(cvsfiles_3, file_Exp3)


# AllDataMerge(c(file_Exp1, file_Exp2, file_Exp3))
# 
# summariseSubInfoAll(c(file_Exp1, file_Exp2, file_Exp3))

AllDataMerge(c(file_Exp1, file_Exp2))

summariseSubInfoAll(c(file_Exp1, file_Exp2))