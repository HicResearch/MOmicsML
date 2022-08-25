#############################################
# File:   extractConfusionMatrix.r
#
# 
# Arguments: output folder and scenarioNumber
#
# Outputs: extracted confusion matrix
#
# Usage: This script extracts confusion matrix from the results of a scenario run.
#
#	Author: PS Reel

extractConfusionMatrix <- function(outDirfinal,scenarioNumber) {
  load(paste(outDirfinal,'Individual_Result_ScenarioNumber_',as.numeric(scenarioNumber),'_output.RData', sep = ''))
  if (exists('output')) {
    t(output$table)
  } else {
    stop('Confusion Matrix not found in the loaded file')
  }
  
}
