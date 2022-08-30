#############################################
# File:   data_preparation.r
#
# 
# Arguments: multi-omics data file in csv format
#
# Outputs: multi-omics data files for training and testing
#
# Usage: This script prepares the multi-omics dataset with disease-wise multi-feature outlier detection and 
# splits the dataset into training and testing sets.
#
#	Author: PS Reel


data_preparation <- function() {

  #load multi-omics file with simulated data
  data <- read.csv('demo_data/Multi-omics/momics_data.csv',fileEncoding="UTF-8-BOM")
  
  c_type <- c('PA','PPGL','CS','PHT')
  
  #find outliers and cap them 
  for (j in 4:(ncol(data)-1)){
    for (i in 1:length(unique(data$Condition))){
      qnt <- quantile(data[which(data$Condition %in% c_type[i]),c(j)], probs=c(.25, .75), na.rm = T)
      caps <- quantile(data[which(data$Condition %in% c_type[i]),c(j)], probs=c(.05, .95), na.rm = T)
      H <- 3 * IQR(data[which(data$Condition %in% c_type[i]),c(j)], na.rm = T)
      data[which(data$Condition %in% c_type[i]),c(j)][data[which(data$Condition %in% c_type[i]),c(j)] < (qnt[1] - H)] <- caps[1]
      data[which(data$Condition %in% c_type[i]),c(j)][data[which(data$Condition %in% c_type[i]),c(j)] > (qnt[2] + H)] <- caps[2]
    }
  }
  #save multi-omics data file without outliers
  write.csv(data, file = 'demo_data/Multi-omics/momics_data_without_outliers.csv', row.names=FALSE)
  
  
  #Split train and test datasets
  library(caTools)
  #random seed for reproducibility
  set.seed(3)
  dataIn <- read.csv('demo_data/Multi-omics/momics_data.csv', stringsAsFactors = TRUE,fileEncoding="UTF-8-BOM")
  
  source('code/stratified.r')
  a<-stratified(dataIn, c("Condition", "Gender"), 0.2)
  dataIn$set <- NA
  dataIn$set[which(dataIn$ID %in% a$ID)] <- "test"
  dataIn$set[is.na(dataIn$set)] <- "train"
  dataIn$set <- factor(dataIn$set, levels = c("train","test"))
  
  write.csv(dataIn[dataIn$set=="train",which(colnames(dataIn)!="set")],"demo_data/Multi-omics/Training/Traindata.csv", row.names = FALSE)
  write.csv(dataIn[dataIn$set=="test",which(colnames(dataIn)!="set")],"demo_data/Multi-omics/Testing/Testdata.csv", row.names = FALSE)

  
  set.seed(3)
  dataIn <- read.csv('demo_data/Multi-omics/momics_data_without_outliers.csv', stringsAsFactors = TRUE,fileEncoding="UTF-8-BOM")
  
  source('code/stratified.r')
  a<-stratified(dataIn, c("Condition", "Gender"), 0.2)
  dataIn$set <- NA
  dataIn$set[which(dataIn$ID %in% a$ID)] <- "test"
  dataIn$set[is.na(dataIn$set)] <- "train"
  dataIn$set <- factor(dataIn$set, levels = c("train","test"))
  
  write.csv(dataIn[dataIn$set=="train",which(colnames(dataIn)!="set")],"demo_data/Multi-omics/Training/Traindata_without_outliers.csv", row.names = FALSE)
  write.csv(dataIn[dataIn$set=="test",which(colnames(dataIn)!="set")],"demo_data/Multi-omics/Testing/Testdata_without_outliers.csv", row.names = FALSE)
}

