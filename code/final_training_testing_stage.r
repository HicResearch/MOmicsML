#############################################
# File:   final_training_testing_stage.r
#
# 
# Arguments: training and testing datasets with ML parameter settings
#
# Outputs: variable
#
# Usage: This script runs the final training testing stage of the biomarker discovery pipeline of the multi-omics analysis.
#
#	Author: PS Reel


final_training_testing_stage <- function() {
  
  
  setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')
  source('classif_Metrics.r')
  
  
  library(RWeka)
  library(caret)
  library(pROC)
     
  #select cut off for feature freq
  n <- 50
  
  datatype <- 'multiomic'
  
  # Training using real training dataset (no SMOTE)
  
  ### Gather top features list
  feature_list <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_1/Set_A/momics/data_for_Venn_diagram_with_cutoff_0.csv")
  Traindata_orig <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/demo_data/Multi-omics/Training/Traindata_without_outliers.csv", stringsAsFactors = TRUE)
  Testdata_orig <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/demo_data/Multi-omics/Testing/Testdata_without_outliers.csv", stringsAsFactors = TRUE)
  
    
  setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Final_Training_Testing_Stage/without_SMOTE/')
  
  ##### PA PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'PA_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PA_PHT_F[(!is.na(feature_list$PA_PHT_C)&(feature_list$PA_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "PA")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "PA")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best classifiers with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
        optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
        classifier_name <- 'SL'
    } else if (i==2) {
        optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
        classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    }  
    
    # chosen randomly to ensure reproducibility  
    set.seed(510000) #510000
    
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PA")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PA")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  
  ##### PPGL PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'PPGL_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PPGL_PHT_F[(!is.na(feature_list$PPGL_PHT_C)&(feature_list$PPGL_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "PPGL")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "PPGL")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PPGL")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PPGL")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  ##### CS PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'CS_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$CS_PHT_F[(!is.na(feature_list$CS_PHT_C)&(feature_list$CS_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "CS")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "CS")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"CS")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"CS")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  ##### ALL ALL ####
  # Select features from final list and subset in Train dataset
  disease <- 'ALL_ALL'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PPGL_PA_CS_PHT_F[(!is.na(feature_list$PPGL_PA_CS_PHT_C)&(feature_list$PPGL_PA_CS_PHT_C >= n))]
  Tdata <- train_data[,c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[,c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PHT") # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PHT")  # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  ##### EHT PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'EHT_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$EHT_PHT_F[(!is.na(feature_list$EHT_PHT_C)&(feature_list$EHT_PHT_C >= n))]
  Tdata <- train_data[,c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- as.character(Tdata$Condition)
  Tdata$Condition[Tdata$Condition =="PA" | Tdata$Condition =="PPGL" | Tdata$Condition =="CS"] <- "EHT"
  Tdata$Condition <- factor(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[,c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- as.character(Vdata$Condition)
  Vdata$Condition[Vdata$Condition =="PA" | Vdata$Condition =="PPGL" | Vdata$Condition =="CS"] <- "EHT"
  Vdata$Condition <- factor(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"EHT") # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"EHT")  # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))

  
  # Training using real & synthetic training dataset (with SMOTE)
  setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/code/')
  
  source('findCommonfeatures.r')
  source('generate_Syndata.r')
  
  feature_out <- findCommonfeatures('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_1/Set_A/momics', 
                                    '//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Final_Training_Testing_Stage/with_SMOTE', 
                                    n)
  generate_Syndata('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/demo_data/Multi-omics/Training/Traindata_without_outliers.csv',
                   '//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Final_Training_Testing_Stage/with_SMOTE', 
                   feature_out[[2]])
  
  feature_list <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Feature_Selection_Stage/without_outliers/Scenarios_subsets/Scenario_1/Set_A/momics/data_for_Venn_diagram_with_cutoff_0.csv")
  Traindata_orig <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Final_Training_Testing_Stage/with_SMOTE/SMOTEdata.csv", stringsAsFactors = TRUE)
  Testdata_orig <- read.csv("//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/demo_data/Multi-omics/Testing/Testdata_without_outliers.csv", stringsAsFactors = TRUE)
  
  setwd('//ensat-vm/shared/Results/5omics_new/Complete_set/code_copy_for_github_release/output/Final_Training_Testing_Stage/with_SMOTE/')
  
  
  ##### PA PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'PA_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PA_PHT_F[(!is.na(feature_list$PA_PHT_C)&(feature_list$PA_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "PA")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "PA")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    }  
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PA")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PA")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  
  ##### PPGL PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'PPGL_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PPGL_PHT_F[(!is.na(feature_list$PPGL_PHT_C)&(feature_list$PPGL_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "PPGL")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "PPGL")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PPGL")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PPGL")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  ##### CS PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'CS_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$CS_PHT_F[(!is.na(feature_list$CS_PHT_C)&(feature_list$CS_PHT_C >= n))]
  Tdata <- train_data[(train_data$Condition== "CS")|(train_data$Condition== "PHT"),c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[(test_data$Condition== "CS")|(test_data$Condition== "PHT"),c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"CS")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"CS")
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  ##### ALL ALL ####
  # Select features from final list and subset in Train dataset
  disease <- 'ALL_ALL'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$PPGL_PA_CS_PHT_F[(!is.na(feature_list$PPGL_PA_CS_PHT_C)&(feature_list$PPGL_PA_CS_PHT_C >= n))]
  Tdata <- train_data[,c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- droplevels(Tdata$Condition)
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[,c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- droplevels(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"PHT") # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"PHT")  # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  ##### EHT PHT ####
  # Select features from final list and subset in Train dataset
  disease <- 'EHT_PHT'
  train_data <-Traindata_orig
  train_data$Condition <- relevel(train_data$Condition,'PPGL') # to fix the order of factors
  features <- feature_list$EHT_PHT_F[(!is.na(feature_list$EHT_PHT_C)&(feature_list$EHT_PHT_C >= n))]
  Tdata <- train_data[,c(which(colnames(train_data) %in%features),which(colnames(train_data)=="Condition")),drop=FALSE]
  Tdata$Condition <- as.character(Tdata$Condition)
  Tdata$Condition[Tdata$Condition =="PA" | Tdata$Condition =="PPGL" | Tdata$Condition =="CS"] <- "EHT"
  Tdata$Condition <- factor(Tdata$Condition)
  set.seed(915) 
  Tdata <- downSample(Tdata, Tdata$Condition)
  Tdata$Class <- NULL
  
  
  test_data <- Testdata_orig
  test_data$Condition <- relevel(test_data$Condition,'PPGL') # to fix the order of factors
  Vdata <- test_data[,c(which(colnames(test_data) %in%features),which(colnames(test_data)=="Condition")),drop=FALSE]
  Vdata$Condition <- as.character(Vdata$Condition)
  Vdata$Condition[Vdata$Condition =="PA" | Vdata$Condition =="PPGL" | Vdata$Condition =="CS"] <- "EHT"
  Vdata$Condition <- factor(Vdata$Condition)
  
  # Train using best one classifier with Training data
  Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
  Train_metrics <- data.frame(matrix(NA, nrow = 1, ncol = 17))
  
  colnames(Train_metrics) <- c('Datatype','Disease','Classifier','feature_freq','feature_count','Accuracy', 'CorrectlyClassified', 'IncorrectlyClassified' ,'TotalInstances', 'Balanced_Accuracy', 'Sensitivity', 'Specificity', 'AUCScore',  'F1', 'Kappa','Precision', 'Recall')
  Test_metrics  <- Train_metrics
  for (i in 1:3){
    if (i==1){
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.functions.SimpleLogistic')))
      classifier_name <- 'SL'
    } else if (i==2) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.trees.RandomForest')))
      classifier_name <- 'RF'
    } else if (i==3) {
      optns <- Weka_control(E = list(c('weka.attributeSelection.InfoGainAttributeEval')), S = list(c('weka.attributeSelection.Ranker')),  W = list(c('weka.classifiers.meta.LogitBoost')))
      classifier_name <- 'LB'
    } 
    
    set.seed(510000) #510000
    trained_Weka_model <- Weka_model(Tdata$Condition ~ ., Tdata,control = optns)
    #Predict classifier with Training data
    pred_out_class <- predict(trained_Weka_model, newdata = Tdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Tdata, type = 'probability')
    metricsT <- classif_Metrics(pred_out_class,pred_out_prob,Tdata$Condition,"EHT") # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Train','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Train','_prob','.csv')) 
    capture.output(metricsT$CM,file =paste0(disease,'_',classifier_name,'_Train','_CM','.csv'))
    Train_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsT[2:13])))
    
    #Predict classifier with Test data
    pred_out_class <- predict(trained_Weka_model, newdata = Vdata) # add type = 'probability' for probabilities instead of labels 
    pred_out_prob <- predict(trained_Weka_model, newdata = Vdata, type = 'probability')
    metricsV <- classif_Metrics(pred_out_class,pred_out_prob,Vdata$Condition,"EHT")  # positive class variable does not effect the outcome for multiclass
    write.csv(pred_out_class,paste0(disease,'_',classifier_name,'_Test','_class','.csv'))
    write.csv(pred_out_prob,paste0(disease,'_',classifier_name,'_Test','_prob','.csv')) 
    capture.output(metricsV$CM,file =paste0(disease,'_',classifier_name,'_Test','_CM','.csv'))
    Test_metrics[i,] <- c(datatype,disease,classifier_name,n,length(features),t(unlist(metricsV[2:13])))
  }
  write.csv(Train_metrics,paste0('Results_Train_',disease,'_Train','.csv')) 
  write.csv(Test_metrics,paste0('Results_Test_',disease,'_Test','.csv')) 
  write.csv(features,paste0('Top_Features_',n,'_',disease,'.csv'))
  
  
  
  
}