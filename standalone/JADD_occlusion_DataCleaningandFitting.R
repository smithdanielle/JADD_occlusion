library(ggplot2)
library(MPDiR)
library(plyr)
library(modelfree)
library(Hmisc)
library(afex)
numTrials <- 150
nReps <- 30

extractData<-function(x){
  data <- read.table(fileName<-x, header=FALSE, skip=60, nrows=numTrials) # import file contents  colnames(data)<-c("stim", "order", "disp", "kresp", "testresp", "correct", "time")
  colnames(data)<-c("stim", "order", "disp", "kresp", "testresp", "correct", "time")
  fileNameSplit <- strsplit(basename(fileName), "_")[[1]] # return name of file only
  info.ID <- as.character(fileNameSplit[1]) # get participant id
  info.disparitySign <- as.character(fileNameSplit[2]) # gets the sign of disparity
  info.condition <- as.character(fileNameSplit[3]) # which condition?
  
  # get data for psychometric function
  intensity<-sort(unique(data$disp))
  nYes<-as.vector(table(data$disp, data$testresp)[,1])
  nNo<-as.vector(table(data$disp, data$testresp)[,2])
  responses <- data.frame(intensity, nYes, nNo)
  responses$ID <- info.ID
  responses$disparitySign<-info.disparitySign
  responses$condition<-info.condition
  
  responses[responses$disparitySign == "uncrossed", c("nYes", "nNo")]<-responses[responses$disparitySign == "uncrossed", c("nNo", "nYes")]
  responses$p <- responses$nYes/(responses$nYes + responses$nNo)
  
  # aggregate median reaction time data
  data.response<-subset(data, , select=c(disp,time))
  data.response.outliersRemoved.full<-subset(data.response, time>0.25 & time<(mean(data.response$time) + 2*sd(data.response$time)))
  data.response.agg.rt<-aggregate(data.response.outliersRemoved.full$time, by=list(disp = data.response.outliersRemoved.full$disp), FUN=median)
  data.response.agg.sd<-aggregate(data.response.outliersRemoved.full$time, by=list(disp = data.response.outliersRemoved.full$disp), FUN=sd)
  data.response.agg.rtReciprocal<-aggregate((1/data.response.outliersRemoved.full$time), by=list(disp = data.response.outliersRemoved.full$disp), FUN=median)
  data.response.agg.sdReciprocal<-aggregate((1/data.response.outliersRemoved.full$time), by=list(disp = data.response.outliersRemoved.full$disp), FUN=sd)
  responses$time<-data.response.agg.rt$x
  responses$sd<-data.response.agg.sd$x
  responses$timeReciprocal<-data.response.agg.rtReciprocal$x
  responses$sdReciprocal<-data.response.agg.sdReciprocal$x
  
  # look at ex-gaussian distrbutions
  library(retimes)
  ExGauss<-function(df){
    results<-mexgauss(df$time)
    results
  }
  rt<-ddply(data.response, .(disp), ExGauss)
  
  responses<-cbind(responses,rt[,2:4])
  
  diag<-data.frame(resp = data$testresp, intensity = data$disp)
  diag$resp[diag$resp == 2]<-0
  diag$ID <- info.ID
  diag$disparitySign<-info.disparitySign
  diag$condition<-info.condition
  
  return(list(responses, diag))
}

my.files <- list.files(pattern = ".txt")
outputs <- lapply(my.files, FUN = extractData)
psychoData <- ldply(outputs, "data.frame")
diagData<-psychoData[,c(15,1,4:6)]
reactionData<-psychoData[,c(1,4:6,8:14)]
colnames(diagData)<-c("resp", "intensity", "ID", "disparitySign", "condition")
colnames(reactionData)<-c("intensity", "ID", "disparitySign", "condition", "rt", "sd", "rtReciprocal", "sdReciprocal", "mu", "sigma", "tau")
psychoData<-psychoData[,1:7]

# Kludge because my code duplicates the results 30 times
psychoData<-psychoData[!duplicated(psychoData),]
reactionData<-reactionData[!duplicated(reactionData),]

# Psychometric function
Fitting<-function(df){
  model<-glm(p ~ intensity, weights = rep(nReps, length(df$intensity)),family = binomial(link = probit), data = df)
  xseq<-seq(min(df$intensity), max(df$intensity), len = 1000)
  yseq<-predict(model, data.frame(intensity = xseq), type = "response")
  
  ID<-unique(df$ID)
  disparitySign<-unique(df$disparitySign)
  conditon<-unique(df$condition)
  
  threshold<-threshold_slope(yseq, xseq, 0.8325)$x_th
  slope<-threshold_slope(yseq, xseq, 0.8325)$slope
  
  data.frame(xseq,yseq,ID,threshold, slope)
}

curves<-ddply(psychoData, .(disparitySign, condition, ID), Fitting)

psychometric.parameters<-aggregate(list(curves$threshold, curves$slope), by=list(curves$ID, curves$disparitySign, curves$condition), 
                      FUN=mean, na.rm=TRUE)
colnames(psychometric.parameters)<-c("ID", "disparitySign", "condition", "threshold", "slope")

# Cleanup for people with backwards psychometric functions
psychometric.parameters$slope<-abs(psychometric.parameters$slope)
psychometric.parameters$threshold<-abs(psychometric.parameters$threshold)

# Assign groups
psychometric.parameters$group<-NA
psychometric.parameters[grep("^A", psychometric.parameters$ID), ]$group<-"asd"
psychometric.parameters[grep("^T", psychometric.parameters$ID), ]$group<-"td"
reactionData$group<-NA
reactionData[grep("^A", reactionData$ID), ]$group<-"asd"
reactionData[grep("^T", reactionData$ID), ]$group<-"td"

# Get stats on disparity levels
# Participants taken out of dataset
occlusion2.remove <- c("A05", "A07", "A19", "A22", "A26", "T16", "T17", "T18", "T27")
dispLevelSummaryData <- psychoData[psychoData$ID %notin% occlusion2.remove, c(1,4,5,6,7)]
summarise(subset(dispLevelSummaryData, intensity != 0), 
          mean = mean (abs(intensity)), median = median(abs(intensity)), 
          mode = as.numeric(names(which.max(table(abs(intensity))))), 
          sd = sd(abs(intensity)), min = min(abs(intensity)), max = max(abs(intensity)))