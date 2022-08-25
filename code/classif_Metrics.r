#############################################
# File:   classif_Metrics.r
#
# 
# Arguments: predicted class labels and probabilities 
#
# Outputs: computed performance metrics
#
# Usage: This script computes perfromance metrics using predicted class labels, probabilities and actual labels.
#
#	Author: PS Reel

classif_Metrics <- function(pred_class,pred_prob,act_class,positive_class)
{
  output <- confusionMatrix(pred_class,act_class, mode='everything', positive = positive_class) 
  ni <- sum(output$table)# number of instances
  nc <- nrow(output$table) # number of classes
  diag <- diag(output$table) # number of correctly classified instances per class 
  
  CorrectlyClassified <- sum(diag)
  IncorrectlyClassified <- ni-CorrectlyClassified
  TotalInstances <- ni
  
  Accuracy <- as.numeric(output$overall['Accuracy'])
  Kappa <- as.numeric(output$overall['Kappa'])
  
  if (nc > 2) { # find averages of class wise metrics in case of multi-class classification
    
    output$byClass <- rbind((output$byClass), colMeans(output$byClass))
    rownames(output$byClass)[nc+1] <- 'Mean Value'
    
    Balanced_Accuracy <- as.numeric(output$byClass['Mean Value','Balanced Accuracy'])
    F1 <- as.numeric(output$byClass['Mean Value','F1'])
    Sensitivity <- as.numeric(output$byClass['Mean Value','Sensitivity'])
    Specificity <- as.numeric(output$byClass['Mean Value','Specificity'])
    Precision <- as.numeric(output$byClass['Mean Value','Precision'])
    Recall <- as.numeric(output$byClass['Mean Value','Recall'])
    
  } else {
    
    Balanced_Accuracy <- as.numeric(output$byClass['Balanced Accuracy'])
    F1 <- as.numeric(output$byClass['F1'])
    Sensitivity <- as.numeric(output$byClass['Sensitivity'])
    Specificity <- as.numeric(output$byClass['Specificity'])
    Precision <- as.numeric(output$byClass['Precision'])
    Recall <- as.numeric(output$byClass['Recall'])
  }
  
  AUC <- auc(multiclass.roc(act_class,pred_prob))[1]
  
  if (is.na(Precision)) warning(paste0('No TP or FP in iteration: ',scenarioIn$scenarioNumber, ', hence Precision set to 0'))
  #https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
  Precision[is.na(Precision)] <- 0 # extreme cases when none classified correctly
  
  if (is.na(Recall)) warning(paste0('No TP or FN in iteration: ',scenarioIn$scenarioNumber, ', hence Recall set to 0'))
  #https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
  Recall[is.na(Recall)] <- 0 # extreme cases when none classified correctly
  
  if (is.na(F1)) warning(paste0('No TP or FN in iteration: ',scenarioIn$scenarioNumber, ', hence F1 set to 0'))
  F1[is.na(F1)] <- 0 # extreme cases when none classified correctly
  
  #classWiseMetrics <- data.frame(Balanced_Accuracy, Kappa, F1, Sensitivity, Specificity, Precision,Recall)
  
  
  return(list('CM' = output$table, 'Accuracy'= Accuracy, 'CorrectlyClassified' = CorrectlyClassified,'IncorrectlyClassified' = IncorrectlyClassified,
              'TotalInstances' = TotalInstances, 'Balanced_Accuracy' = Balanced_Accuracy, 'Sensitivity' = Sensitivity, 'Specificity' = Specificity,
              'AUCScore' = AUC, 'F1' = F1, 'Kappa' = Kappa,'Precision' = Precision,'Recall' = Recall))  
  
}