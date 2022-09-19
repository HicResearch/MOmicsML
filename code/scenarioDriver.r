#############################################
# File:   scenarioDriver.r
#
# 
# Arguments: various input parameters (see below):
#  fileName_Entries: dataset file path
#  output_Directory_Entries: file path to save results
#  omicsData_Entries: Omics type
#  clinicalPhenotype_Entries: clinical Phenotypes to include
#  conditionToClassify_Entries: disease combinations to include 
#  validationType_Entries: validation type to perform
#  classifierMethod_Entries: type of classifer and feature selection method to use
#  randomSeed_Entries: number of random repeats
#  feature_preselection_flag: use preselected features
#  feature_preselectionpath: file path of preseelcted features file
#  feature_preselection_freqcutoff: freq cutoff for features
#  gensyn_flag: allow syn generation for training set 
#
#
# Outputs: performance metrics for scenario runs
#
# Usage: This script drives the scenarioRunner function which runs each random repeat.
#
#	Author: PS Reel

scenarioDriver <- function(input_parameters) {
  
  # Increase the JVM heap size in rJava's options support:
  # options(java.parameters = "-Xmx2048m") # if you want 2GB RAM
  options(java.parameters = "-Xmx2048m")
  
  library(caret)
  library(caTools)
  library(RWeka)
  library(rJava)
  library(plyr)
  library(stringr)
  library(Boruta)
  library(pROC)
  
  source("code/findCommonfeatures.r")
  source("code/generate_Syndata.r")
  source("code/scenarioRunner.r")
  source("code/scenarioGenerator.r")
  source("code/extractFeatures.r")
  source("code/generate_CM.r")
  source("code/scenarioResultAnalysis_Summary.r")
  
  
  
  allScenarios <- scenarioGenerator(input_parameters)
  
  zz <- file(paste0(input_parameters$output_Directory_Entries, "/Console_output_file.txt"), open="wt")
  #sink(paste0(input_parameters$output_Directory_Entries, "/Console_output_file.txt"))
  #sink(zz ,type = "output",split = F)
  #sink(zz, append = T, type = "output", split = T) 
  
  
  if (input_parameters$feature_preselection_flag == 1) {
    print('Feature preselection ENABLED')
    feature_out <- findCommonfeatures(input_parameters$feature_preselectionpath, input_parameters$output_Directory_Entries, input_parameters$feature_preselection_freqcutoff)
  }
  
  if ((input_parameters$gensyn_flag == 1) & (input_parameters$feature_preselection_flag = 0)) {
    stop('Do not generate SMOTE dataset on complete dataset. Please preselect features using a freq file')
  }
  
  if (input_parameters$gensyn_flag == 1) {
    print('SMOTE based Synthetic data generation ENABLED')
    generate_Syndata(input_parameters$fileName_Entries,input_parameters$output_Directory_Entries, feature_out[[2]])
  }
  
 
  
  write.csv(allScenarios,paste0(input_parameters$output_Directory_Entries, "/input_scenarios_summary.csv"),row.names = FALSE)
  print('Input Classification Scenarios saved')
  numOfScenarios <- nrow(allScenarios)
  
  allScenariosSummaryResult <- allScenarios
  allScenariosSummaryResult$Accuracy <- NA
  allScenariosSummaryResult$CorrectlyClassified <- NA
  allScenariosSummaryResult$IncorrectlyClassified <- NA
  allScenariosSummaryResult$TotalInstances <- NA
  allScenariosSummaryResult$Balanced_Accuracy <- NA
  allScenariosSummaryResult$Kappa <- NA
  allScenariosSummaryResult$F1 <- NA
  allScenariosSummaryResult$Sensitivity <- NA
  allScenariosSummaryResult$Specificity <- NA
  allScenariosSummaryResult$Precision <- NA
  allScenariosSummaryResult$Recall <- NA
  allScenariosSummaryResult$ComputeTime <- NA
  allScenariosSummaryResult$AUCScore <- NA
  
  
     
  for (i in 1:nrow(allScenarios)) { # 1:nrow(allScenarios)
    temp_results <- scenarioRunner(allScenarios[i,])
    allScenariosSummaryResult$Accuracy[i] <- temp_results$Accuracy
    allScenariosSummaryResult$CorrectlyClassified[i] <- temp_results$CorrectlyClassified
    allScenariosSummaryResult$IncorrectlyClassified[i] <- temp_results$IncorrectlyClassified
    allScenariosSummaryResult$TotalInstances[i] <- temp_results$TotalInstances
    allScenariosSummaryResult$Balanced_Accuracy[i] <- temp_results$Balanced_Accuracy
    allScenariosSummaryResult$Kappa[i] <- temp_results$Kappa
    allScenariosSummaryResult$F1[i] <- temp_results$F1
    allScenariosSummaryResult$Sensitivity[i] <- temp_results$Sensitivity
    allScenariosSummaryResult$Specificity[i] <- temp_results$Specificity
    allScenariosSummaryResult$Precision[i] <- temp_results$Precision
    allScenariosSummaryResult$Recall[i] <- temp_results$Recall
    allScenariosSummaryResult$ComputeTime[i] <- temp_results$ComputeTime
    allScenariosSummaryResult$AUCScore[i] <- temp_results$AUCScore
    
    print(paste("Successfully processed",i, "out of",numOfScenarios,".",  sep=" "))
  }
  
  
  
  write.csv(allScenariosSummaryResult, file = paste(allScenarios$outDir[1],"allScenariosSummaryResult.csv", sep=""), row.names = FALSE)
  
  
  #extract selected features from all scenarios and append to file
  extractFeatures(allScenarios$outDir[1])
  
  generate_CM(allScenarios$outDir[1])
  
  #generate perfromance and selected features plots
  scenarioResultAnalysis_Summary(allScenarios$outDir[1])
  
  
  print('Scenario completed successfully')
 
  
}