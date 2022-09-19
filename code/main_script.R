#############################################
# File:   main_script.r
#
# 
# Arguments: variable
#
# Outputs: variable
#
# Usage: This script runs the complete biomarker discovery pipeline of the multi-omics analysis.
#
#	Author: PS Reel

#setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')

# Data Preparation
source('code/data_preparation.r')
data_preparation()

# Preprocessing Stage
source('code/preprocessing_stage.r')
preprocessing_stage()

# Feature Selection Stage
source('feature_selection_stage.r')
feature_selection_stage()

# Final Training Testing Stage
source('final_training_testing_stage.r')
final_training_testing_stage()
