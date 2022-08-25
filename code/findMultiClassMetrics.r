#############################################
# File:   findMultiClassMetrics.r
#
# 
# Arguments: multiclass confusion matrix 
#
# Outputs: computed performance metrics
#
# Usage: This script computes perfromance metrics from multiclass confusion matrix.
#
#	Author: PS Reel

findMultiClassMetrics <- function(CM) {

#find different metrics from multiclass confusion matrix
# such as accuracy, precision, specificity, sensitivity and F-score
# https://www.researchgate.net/post/Can_someone_help_me_to_calculate_accuracy_sensitivity_of_a_66_confusion_matrix

  tpn <- sum(diag(CM))
  tp <- t(diag(CM))
  
  tn <- fn <- fp <- NA
  acc <- prec <- sen <- spec <- fscor <- NA
  
  for (i in 1:dim(CM)[1]){
    tn[i] <- tpn - tp[i]
    fn[i] <- sum(CM[i,]) - CM[i,i]
    fp[i] <- sum(CM[,i]) - CM[i,i]
  }
   
  for (i in 1:dim(CM)[1]){
    acc[i] <- ((tp[i] + tn[i]) / (tp[i] + fp[i] + fn[i] + tn[i])) * 100
    
    prec[i] <- (tp[i] / (tp[i] + fp[i])) * 100
    
    
    sen[i] <- (tp[i] / (tp[i] + fn[i])) * 100
    
    
    spec[i] <- (tn[i] / (fp[i] + tn[i])) * 100
    
    
    fscor[i] <- (sen[i] * prec[i] * 2) / (sen[i] + prec[i])
  } 
  
  prec[(tp==0) & (fp==0) & (fn==0)] <- 100
  prec[(tp==0) & ((fp!=0) | (fn!=0))] <- 0
  sen[(tp==0) & (fp==0) & (fn==0)] <- 100
  sen[(tp==0) & ((fp!=0) | (fn!=0))] <- 0
  spec[(tn==0) & (fp==0) & (fn==0)] <- 100
  spec[(tn==0) & (fn==0)] <- 0
  
  overall_acc <- (tpn/sum(sum(CM)))*100
  
  W_prec <- sum(prec*rowSums(CM))/sum(rowSums(CM)) 
  W_sen <- sum(sen*rowSums(CM))/sum(rowSums(CM)) 
  W_spec <- sum(spec*rowSums(CM))/sum(rowSums(CM)) 
  W_fscor <- sum(fscor*rowSums(CM))/sum(rowSums(CM)) 

  out_data <- list('W_sen'= W_sen/100,'W_spec'= W_spec/100)
}