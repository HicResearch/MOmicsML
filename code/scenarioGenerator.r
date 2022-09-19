#############################################
# File:   scenarioGenerator.r
#
# 
# Arguments: various input parameters (see below):
#
# Outputs: performance metrics for scenario runs
#
# Usage: This script prepares the instructions needed to run each scenario.
#
#	Author: PS Reel

scenarioGenerator <- function(input_parameters) {

  fileName_Entries <- input_parameters$fileName_Entries
  output_Directory_Entries <- input_parameters$output_Directory_Entries
  omicsData_Entries <- input_parameters$omicsData_Entries
  clinicalPhenotype_Entries <- input_parameters$clinicalPhenotype_Entries
  conditionToClassify_Entries <- input_parameters$conditionToClassify_Entries
  validationType_Entries <- input_parameters$validationType_Entries
  classifierMethod_Entries <- input_parameters$classifierMethod_Entries
  randomSeed_Entries <- input_parameters$randomSeed_Entries
  feature_preselection_flag <- input_parameters$feature_preselection_flag
  gensyn_flag <- input_parameters$gensyn_flag
  
  # Check input file if it exists
  if (file.exists(fileName_Entries)){
    print('The input data file exists.')
  } else {
    stop('The input file: ',fileName_Entries , ' does NOT exist. Please check the file path and run again.')
  }
  
  # check for forward slash at the end
  if (substr(output_Directory_Entries,nchar(output_Directory_Entries),nchar(output_Directory_Entries)) == '/'){
    stop('Remove last / from the output folder path and run again.')
  }
  
  # Check if the output folder exists
  if (!file.exists(output_Directory_Entries)){
    print('The output folder does not exist. Creating now...')
    dir.create(file.path(output_Directory_Entries), recursive = TRUE)
    output_Directory_Entries <- paste0(output_Directory_Entries,'/')
  } else {
    stop(paste0('The output folder: ',output_Directory_Entries , ' already EXISTS. Please check the output folder path and run again.'))
  }
  
  #write the scenario config
  sink(paste0(output_Directory_Entries,"config_file.txt"))
  print(input_parameters)
  sink()
  
  totalScenarios <- expand.grid(fileName_Entries, output_Directory_Entries,omicsData_Entries, clinicalPhenotype_Entries, conditionToClassify_Entries, validationType_Entries, 
                                classifierMethod_Entries, randomSeed_Entries, feature_preselection_flag, gensyn_flag)
  totalScenarios <- cbind(scenarioNumber <- c(1:nrow(totalScenarios)), totalScenarios)
  colnames(totalScenarios) <- c("scenarioNumber" ,"fileName", "outDir","omicsData", "clinicalPhenotype", "conditionToClassify", "validationType", 
                                "classifierMethod", "randomSeed", "feature_preselection_flag","gensyn_flag")
  totalScenarios
}