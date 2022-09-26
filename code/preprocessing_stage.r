#############################################
# File:   preprocessing_stage.r
#
# 
# Arguments: training dataset (with and without outliers) with ML parameter settings
#
# Outputs: variable
#
# Usage: This script runs the preprocessing stage of the biomarker discovery pipeline of the multi-omics analysis.
#
#	Author: PS Reel

preprocessing_stage <- function() {

  source("code/scenarioDriver.r")
  
  #with outliers All multi-omic features
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata.csv"),
    output_Directory_Entries = c("output/Preprocessing_Stage/with_outliers/All_features/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"),
    conditionToClassify_Entries = c("PPGL_PA_CS_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Infogain_RankerAll_RF",  "RWeka_ASC_Infogain_RankerAll_J48", "RWeka_ASC_Infogain_RankerAll_LMT", "RWeka_ASC_Infogain_RankerAll_SL",
                                 "RWeka_ASC_Infogain_RankerAll_SMO", "RWeka_ASC_Infogain_RankerAll_NB",  "RWeka_ASC_Infogain_RankerAll_IBk",
                                 "RWeka_ASC_Infogain_RankerAll_LB"),
    randomSeed_Entries = c(1:100),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  #without outliers All multi-omic features
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Preprocessing_Stage/without_outliers/All_features/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"),
    conditionToClassify_Entries = c("PPGL_PA_CS_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Infogain_RankerAll_RF",  "RWeka_ASC_Infogain_RankerAll_J48", "RWeka_ASC_Infogain_RankerAll_LMT", "RWeka_ASC_Infogain_RankerAll_SL",
                                 "RWeka_ASC_Infogain_RankerAll_SMO", "RWeka_ASC_Infogain_RankerAll_NB",  "RWeka_ASC_Infogain_RankerAll_IBk",
                                 "RWeka_ASC_Infogain_RankerAll_LB"),
    randomSeed_Entries = c(1:100),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
}

# Compare the performance metrics (Accuracy, Sensitivity and Specificity) for the above to select dataset type and best 3 classifiers.