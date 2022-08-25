#############################################
# File:   findCommonfeatures.r
#
# 
# Arguments: input, output folder path and feature freq cut off
#
# Outputs: list of common features across all disease combinations
#
# Usage: This script finds the common features across different disease combinations from scenario runs.
#
#	Author: PS Reel

findCommonfeatures <- function(data_file_inpath, data_file_outpath, cutoff) {
# input is feature list and  cut off
# output is common and individual list of features
# this will be used for SMOTE data generation
  library(tidyr)
  library(memisc)
  
  features <-  read.csv(paste0(data_file_inpath,'/data_for_Venn_diagram_with_cutoff_0.csv')) # give path for now
  
  
  for (i in  seq( 2, ncol(features), 2)) {
    temp1<- which(features[[i]] < cutoff)   
    features[temp1,i-1] <- ''
    features[temp1,i] <- NA
  }                         
  
  
  feature_list <- features[-seq(2,ncol(features),2)]
  
  
  commonfeature_list <- gather(feature_list)      
  commonfeature_list <- commonfeature_list[-c(1)]

  commonfeature_list[commonfeature_list==""] <- NA
  commonfeature_list <- na.omit(commonfeature_list)
  commonfeature_list <- unique(commonfeature_list)
  
  
  commonfeature_list$value <- sub("Gender \\{Female,Male\\}", "Gender", commonfeature_list$value)
  commonfeature_list$value <- sub("O1_VisibleHaemolysis \\{No,Yes\\}", "VisibleHaemolysis", commonfeature_list$value)
  
  
  commonfeature_list <- rbind(commonfeature_list, 'Age','Gender','Condition') 
  commonfeature_list <- unique(commonfeature_list)
  
  commonfeature_list <- as.character(commonfeature_list[,1])
  
  out<- list(features,commonfeature_list)
  write.csv(features, file = paste0(data_file_outpath,"/individual_featurelist.csv"), row.names = FALSE, na = '')
  print('Common features computed for SMOTE generation successfully')
  out
} 






#   ## Compare and drop unwanted features
#   list <- as.character(final_feature_list[,1])
#   var.out.bool <- names(dat) %in% list
#   dat <- dat[,var.out.bool]
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
# feature_list <- read.csv('\\\\ensat-vm/shared/Results for monoomic manuscripts/miRNA/paper pipeline/Step 1 Outliers/without outliers/Boruta/data_for_Venn_diagram_with_cutoff_0.csv')
# 
# feature_cols <- c(paste0(as.character(scenarioIn$conditionToClassify),'_F'), paste0(as.character(scenarioIn$conditionToClassify),'_C'))
# selected_features <- feature_list[ ,which(names(feature_list) %in% feature_cols)]
# 
# # make in individual based on each dissease combination, if needed
# if (as.character(scenarioIn$conditionToClassify) == "CS_PHT"){
#   feature_freq_cutoff <- 50
# } else if (as.character(scenarioIn$conditionToClassify) == "EHT_PHT"){
#   feature_freq_cutoff <- 50
# } else if (as.character(scenarioIn$conditionToClassify) == "PA_PHT"){
#   feature_freq_cutoff <- 50
# } else if (as.character(scenarioIn$conditionToClassify) == "PPGL_PA_CS_PHT"){
#   feature_freq_cutoff <- 50
# } else if (as.character(scenarioIn$conditionToClassify) == "PPGL_PHT"){
#   feature_freq_cutoff <- 50
# }
# 
# pre_feature_list <- as.character(selected_features[which(selected_features[2] >= feature_freq_cutoff), 1])
# print(paste0("Features selected for ", as.character(scenarioIn$conditionToClassify), " which appear ", feature_freq_cutoff ," times or greater"))
# #check for nominal variable such as Age and Visible Heomoysis and fix them
# if (length(grep(c("Gender"), (pre_feature_list))) != 0){
#   pre_feature_list[grep(c("Gender"), (pre_feature_list))] <- "Gender"
# }
# #do similar for hemolysis
# 
# feature_exclude <- c(grep(c("Age"), (pre_feature_list)), grep(c("Gender"), (pre_feature_list)))
# if (omic_prefix == "Multiomic"){
#   stop("Fix the pre selected feature selection part")
#   
# } else {
#   if (length(feature_exclude) != 0){
#     pre_feature_list[-feature_exclude] <- paste0(omic_prefix, pre_feature_list[-feature_exclude]) 
#   } else {
#     pre_feature_list <- paste0(omic_prefix, pre_feature_list) 
#   }
#   
# }
# 
# #which(pre_feature_list %in% colnames(dataIn))
# dataIn <- dataIn[,unlist(lapply(c(pre_feature_list, "Condition"), function(x) grep(paste0("\\b",x,"\\b"), colnames(dataIn))))]
# #dim(dataIn[ ,which(names(dataIn) %in% c(as.character(pre_feature_list), "Condition"))])
# } else {
#   print("feature_preselection_flag is DISABLED. NO feature preselection in employed")
# }
