#############################################
# File:   scenario_subset_generator.r
#
# 
# Arguments: training data path
#
# Outputs: training data subsets
#
# Usage: This script splits the training data into subsets based on Sex (Set C: Males & D: Females) and Age (Set E: Age>=50 and F: Age < 50 years).
#
#	Author: PS Reel

scenario_subset_generator <- function(training_data_path) {
  
  setwd(training_data_path)
  file_name="Traindata_without_outliers"
  
  all_sample<- read.csv(paste0(file_name,".csv"))
  #Age subsets
  for (i in 1:2){
    if (i==1){#male scenario
      age_type='A50'
      subset <- all_sample[which(all_sample$Age >= 50),]
    } else if (i==2){#female scenario
      age_type='B50'
      subset <- all_sample[which(all_sample$Age < 50),]
    }
    
    subset <- droplevels(subset)
    write.csv(subset, file = paste0("Scenario_subsets/Age_subset/",age_type,"/",file_name,".csv"), row.names=FALSE)
  }
  
  #Sex subsets
  for (i in 1:2){
    if (i==1){#male scenario
      gender_type='Male'
    } else if (i==2){#female scenario
      gender_type='Female'
    }
    
    subset <- all_sample[which(all_sample$Gender %in% gender_type),]
    subset <- droplevels(subset)
    write.csv(subset, file = paste0("Scenario_subsets/Sex_subset/",gender_type,"/",file_name,".csv"), row.names=FALSE)
    
  }
  
}