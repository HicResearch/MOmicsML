#############################################
# File:   generate_Syndata.r
#
# 
# Arguments: input, output file path and list of common features
#
# Outputs: dataframe of syn generated data
#
# Usage: This script generates syn data from real data using SMOTE method.
#
#	Author: PS Reel

generate_Syndata <- function(data_file_inpath, data_file_outpath, commonfeature_list) {
  library(DMwR2)
  library(UBL)
  library(sqldf)
  orig_data<- read.csv(paste0(data_file_inpath), stringsAsFactors=TRUE) # give path for now
  
   
  orig_data_subset <- orig_data[,which(colnames(orig_data) %in% commonfeature_list)]
  orig_data_subset <- droplevels(orig_data_subset)
  
  print(table(orig_data_subset$Condition))
  
    
  orig_data_subset$Age <- as.numeric(orig_data_subset$Age)
  orig_data_subset$Gender <- as.numeric(orig_data_subset$Gender)
  
  if ("O1_VisibleHaemolysis" %in% colnames(orig_data_subset)){
    orig_data_subset$O1_VisibleHaemolysis <- as.numeric(orig_data_subset$O1_VisibleHaemolysis)
    }
  
  #For balanced or extreme data set generation
  SMOTE_data <- SmoteClassif(Condition ~ ., orig_data_subset, C.perc =  list(CS = 3, PA = 1, PHT = 1, PPGL =1.5), k = 5, repl = FALSE, dist = "Euclidean", p = 2)
  print(table(SMOTE_data$Condition))
  
  #Mark the original and new rows
  old <- sqldf('SELECT * FROM orig_data_subset INTERSECT SELECT * FROM SMOTE_data')
  SMOTE_new <- sqldf('SELECT * FROM SMOTE_data EXCEPT SELECT * FROM orig_data_subset')
    
  old$datatype <- 1
  SMOTE_new$datatype <- 2
    
  SMOTE_all <- rbind(old,SMOTE_new)
    
  write.csv(SMOTE_all, file = paste0(data_file_outpath,"/SMOTEdata.csv"), row.names = FALSE)
  
  print('SMOTE data generated successfully')
  
  }