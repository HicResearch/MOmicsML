#############################################
# File:   feature_selection_stage.r
#
# 
# Arguments: training dataset with ML parameter settings
#
# Outputs: performance metrics for various scenarios
#
# Usage: This script runs the feature selection stage of the biomarker discovery pipeline of the multi-omics analysis.
#
#	Author: PS Reel

feature_selection_stage <- function() {
  
  #setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')
  source("code/scenarioDriver.r")
  
  #without outliers CFS
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/CFS/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"),
    conditionToClassify_Entries = c("PPGL_PA_CS_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_CFS_Best@F_RF",  "RWeka_ASC_CFS_Best@F_SL", "RWeka_ASC_CFS_Best@F_LB"),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()


  #without outliers Boruta
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Boruta/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"),
    conditionToClassify_Entries = c("PPGL_PA_CS_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF",  "RWeka_ASC_Boruta_SL", "RWeka_ASC_Boruta_LB"),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()

  #without outliers All features
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/All_features/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"),
    conditionToClassify_Entries = c("PPGL_PA_CS_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Infogain_RankerAll_RF", "RWeka_ASC_Infogain_RankerAll_SL", "RWeka_ASC_Infogain_RankerAll_LB"),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  # Compare the performance metrics (Accuracy, Sensitivity and Specificity) for the above to select the best feature selection method.
  
  
  # Following scenarios include the use of age and sex as features and understanding the effect of age and sex segregated subsets on 
  # feature selection in different disease combinations
  
  # Generate data subsets for different scenarios
  #setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')
  source("code/scenario_subset_generator.r")
  scenario_subset_generator('demo_data/Multi-omics/Training/')
  
  #Scenario 1: Including (Set A) vs excluding (Set B) age and sex as features
  #setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')
  # Set A
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_1/Set_A/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("all"), # includes Age and Sex
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","CS_PHT", "PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  # Set B
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_1/Set_B/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("none"), # excludes Age and Sex 
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","CS_PHT", "PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  
  #Scenario 2: Males (Set C) vs Females (Set D)
  
  #Set C
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Scenario_subsets/Sex_subset/Male/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_2/Set_C/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("ageOnly"),
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  #Set D
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Scenario_subsets/Sex_subset/Female/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_2/Set_D/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("ageOnly"),
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","CS_PHT", "PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  #Scenario 3: Older (Set E) vs Younger (Set F)
  
  #Set E
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Scenario_subsets/Age_subset/A50/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_3/Set_E/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("genderOnly"),
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","CS_PHT", "PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  #Set F
  scenarioDriver(list(
    fileName_Entries = c("demo_data/Multi-omics/Training/Scenario_subsets/Age_subset/B50/Traindata_without_outliers.csv"),
    output_Directory_Entries = c("output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_3/Set_F/momics"),
    omicsData_Entries = c('all'),
    clinicalPhenotype_Entries = c("genderOnly"),
    conditionToClassify_Entries = c("PA_PHT","PPGL_PHT","CS_PHT", "PPGL_PA_CS_PHT","EHT_PHT"),
    validationType_Entries = c("DS_80_20"),
    classifierMethod_Entries = c("RWeka_ASC_Boruta_RF","RWeka_ASC_Boruta_SL",'RWeka_ASC_Boruta_LB'),
    randomSeed_Entries = c(1:2),
    feature_preselection_flag = 0,
    feature_preselectionpath = c('null'),
    feature_preselection_freqcutoff = 0,
    gensyn_flag = 0
  ))
  warnings()
  sink(type="output")
  sink(type="message")
  closeAllConnections()
  
  # Compare the performance metrics (Accuracy, Sensitivity and Specificity) for the above scenarios and save the the top selected features.
  
}