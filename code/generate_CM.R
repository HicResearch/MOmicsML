#############################################
# File:   generate_CM.r
#
# 
# Arguments: results output folder path 
#
# Outputs: confusion matrix
#
# Usage: This script generates the confusion matrix from the scenario runs.
#
#	Author: PS Reel

generate_CM <- function(outDirfinal) {
  library(plyr)
  source("extractConfusionMatrix.r")
  source("findMultiClassMetrics.r")
  source("saveConfusionMatrix.r")
  

  scenarioData <- read.csv(paste(outDirfinal,'allScenariosSummaryResult.csv', sep=""),stringsAsFactors = TRUE)
  test <- data.frame(1,1)
  
   #parallel code
   library(future.apply)

   plan('multisession')
   future_lapply(1:nrow(scenarioData), function(i) future({saveConfusionMatrix(outDirfinal,i,extractConfusionMatrix(outDirfinal,i))}),future.seed = TRUE)
    
  
  # code to find the best scenarios for each classification setup and save it in another folder
  
  dir.create(file.path(outDirfinal, 'ConfusionMatrix','best'))
  for(k in 1:length(levels(scenarioData$conditionToClassify))){
    
    scenarioData_subset <- scenarioData[which(scenarioData$conditionToClassify %in% levels(scenarioData$conditionToClassify)[k]),]
    best_scenario <-scenarioData_subset$scenarioNumber[which(scenarioData_subset$Balanced_Accuracy == max(scenarioData_subset$Balanced_Accuracy), arr.ind = TRUE)[1]] # select only the first of more than one
    
    file.copy(file.path(outDirfinal, 'ConfusionMatrix',paste0('CM_',best_scenario,'.png')), file.path(outDirfinal, 'ConfusionMatrix','best', paste0('CM_',best_scenario,'_',levels(scenarioData$conditionToClassify)[k] ,'.png')))
  }
  
  }
