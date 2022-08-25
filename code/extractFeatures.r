#############################################
# File:   extractFeatures.r
#
# 
# Arguments: output results folder
#
# Outputs: extracted features
#
# Usage: This script extracts the features from all scenario runs.
#
#	Author: PS Reel

extractFeatures <- function(outDirfinal) {
  library(plyr)
  scenarioData <- read.csv(paste(outDirfinal,'allScenariosSummaryResult.csv',sep=""))
  
  temp_frame <- data.frame(t(paste(rep("test_",5),c(1:5),sep="")))
  
  colnames(temp_frame)<-paste(rep("Feature_",ncol(temp_frame)),c(1:ncol(temp_frame)),sep="")
  rownames(temp_frame)[1] <- 100
  
  for (i in 1:nrow(scenarioData)){ # 
    
    print(i)
    temp_filename <- paste(outDirfinal,'Individual_Result_ScenarioNumber_',as.numeric(scenarioData$scenarioNumber[i]),'_output.RData', sep = '')
    if (file.exists(temp_filename)) {
    load(temp_filename)
    
    textfree_Weka <- trained_Weka_model$classifier$toString()
    textframe_Weka <- data.frame(strsplit(textfree_Weka, split = '\n'), stringsAsFactors = FALSE)
    colnames(textframe_Weka)[1] <- 'textdata'
    
    
    text_temp <- textframe_Weka[grep("@attribute", textframe_Weka$textdata),1] 
    final_frame <- data.frame(t(text_temp[-grep("@attribute dataTrain", text_temp)])) 
    
    } else {
      final_frame <- data.frame(NA)
    }
    
    colnames(final_frame)<-paste(rep("Feature_",ncol(final_frame)),c(1:ncol(final_frame)),sep="")
    
    
    l <- list(temp_frame, final_frame)
    temp_frame <- rbind.fill(temp_frame,final_frame)
    rm(list=setdiff(ls(), c("i","outDirfinal","scenarioData","temp_frame")))
    gc()
  }
  temp_frame <- temp_frame[-1,]
  rownames(temp_frame)[1:nrow(temp_frame)] <- as.numeric(rownames(temp_frame)[1:nrow(temp_frame)])-1
  
  temp_frame[] <- lapply(temp_frame, gsub, pattern='@attribute ', replacement='')
  temp_frame[] <- lapply(temp_frame, gsub, pattern=' numeric', replacement='')  
  temp_frame[] <- lapply(temp_frame, gsub, pattern=' nominal', replacement='')  
  
  scenarioData <- cbind(scenarioData,temp_frame, row.names=NULL)
  write.csv(scenarioData, file = paste(outDirfinal,'allScenariosSummaryResult.csv', sep=""), row.names = FALSE, na = '')
  
  
  print('Results file is successfully appended with selected features.')
  
 }
