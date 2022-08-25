#############################################
# File:   scenarioRunner.r
#
# 
# Arguments: various input parameters for individual scenario run(see below):
#  fileName_Entries: dataset file path
#  output_Directory_Entries: file path to save results
#  omicsData_Entries: Omics type
#  clinicalPhenotype_Entries: clinical Phenotypes to include
#  conditionToClassify_Entries: disease combinations to include 
#  validationType_Entries: validation type to perform
#  classifierMethod_Entries: type of classifer and feature selection method to use
#  randomSeed_Entries: number of random repeats
#  feature_preselection_flag: use preselected features
#  feature_preselectionpath: file path of preseelcted features file
#  feature_preselection_freqcutoff: freq cutoff for features
#  gensyn_flag: allow syn generation for training set #
#
# Outputs: performance metrics for scenario runs
#
# Usage: This script runs the individual scenario.
#
#	Author: PS Reel

scenarioRunner <- function(scenarioIn) {
  
  if (scenarioIn$gensyn_flag == 1) {
    ensat_ht_data <- data.frame(read.csv(paste0(as.character(scenarioIn$outDir),'/SMOTEdata.csv'),stringsAsFactors = TRUE))
    ensat_ht_data$Condition <- relevel(ensat_ht_data$Condition,'PPGL')
    } else {
    ensat_ht_data <- data.frame(read.csv(as.character(scenarioIn$fileName),stringsAsFactors = TRUE))
    ensat_ht_data$Condition <- relevel(ensat_ht_data$Condition,'PPGL')
    }
  
  if (scenarioIn$feature_preselection_flag == 1) { 
    preselected_features <- read.csv(paste0(as.character(scenarioIn$outDir),"/individual_featurelist.csv"))
    preselected_features <- preselected_features[ , -c(grep("T_C", colnames(preselected_features)))]
  }
   
  if (scenarioIn$feature_preselection_flag == 0) {  
    ### choose omics ####
    switch(as.character(scenarioIn$omicsData), 
           Pmirna = { omic_prefix <- "O1_"
               dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep(omic_prefix, colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           Pmetas = { omic_prefix <- "O2_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O2_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           PSteroids = { omic_prefix <- "O3_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O3_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           USteroids = { omic_prefix <- "O4_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O4_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           SmallMetabolites = { omic_prefix <- "O5_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O5_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           PNMR = { omic_prefix <- "O6_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O6_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           UNMR = { omic_prefix <- "O7_"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O7_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           all ={ omic_prefix <- "Multiomic"
             dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("O1_", colnames(ensat_ht_data)),grep("O2_", colnames(ensat_ht_data)),grep("O3_", colnames(ensat_ht_data)),grep("O4_", colnames(ensat_ht_data)),grep("O5_", colnames(ensat_ht_data)),grep("O6_", colnames(ensat_ht_data)),grep("O7_", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]},
           Methy ={ omic_prefix <- "Methy"
           dataIn <- ensat_ht_data[,c(grep("Gender", colnames(ensat_ht_data)),grep("Age", colnames(ensat_ht_data)),grep("cg", colnames(ensat_ht_data)),which(colnames(ensat_ht_data)=="Condition")),drop=FALSE]}
           ) 
} else {
  dataIn <- ensat_ht_data
  }
    
    
  if (scenarioIn$feature_preselection_flag == 0) {
    print("feature_preselection_flag is DISABLED. NO feature preselection in employed")
    # choose clinical phenotypes ####
    switch(as.character(scenarioIn$clinicalPhenotype), 
           genderOnly = {dataIn <- dataIn[,-c(grep("Age", colnames(dataIn))),drop=FALSE]},
           ageOnly ={dataIn <- dataIn[,-c(grep("Gender", colnames(dataIn))),drop=FALSE]},
           all ={dataIn <- dataIn},
           none ={dataIn <- dataIn[,-c(grep("Age", colnames(dataIn)),grep("Gender", colnames(dataIn))),drop=FALSE]} )
  }
  
  if (scenarioIn$feature_preselection_flag == 1) {
    print('feature_preselection_flag is ENABLED. Hence ONLY pre-selected features used')
    if (scenarioIn$gensyn_flag == 0){
      dataIn <- dataIn[,which(colnames(dataIn) %in% c(as.character(preselected_features[,grep(paste0('^',scenarioIn$conditionToClassify,'_F'),colnames(preselected_features))]) ,'Condition')) ] 
      }
    if (scenarioIn$gensyn_flag == 1){
      dataIn <- dataIn[,which(colnames(dataIn) %in% c(as.character(preselected_features[,grep(paste0('^',scenarioIn$conditionToClassify,'_F'),colnames(preselected_features))]) ,'Condition','datatype')) ] 
      }
    }
          
  
  
    # choose condition to classify 
    switch(as.character(scenarioIn$conditionToClassify), 
           PPGL_PA = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PA'),]},
           PPGL_CS = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='CS'),]},
           PPGL_PHT = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PHT'),]},
           PPGL_HV = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='HV'),]},
           PA_CS = {dataIn <- dataIn[which(dataIn$Condition=='PA'|dataIn$Condition=='CS'),]},
           PA_PHT = {dataIn <- dataIn[which(dataIn$Condition=='PA'|dataIn$Condition=='PHT'),]},
           PA_HV = {dataIn <- dataIn[which(dataIn$Condition=='PA'|dataIn$Condition=='HV'),]},
           CS_PHT = {dataIn <- dataIn[which(dataIn$Condition=='CS'|dataIn$Condition=='PHT'),]},
           CS_HV = {dataIn <- dataIn[which(dataIn$Condition=='CS'|dataIn$Condition=='HV'),]},
           PHT_HV = {dataIn <- dataIn[which(dataIn$Condition=='PHT'|dataIn$Condition=='HV'),]},
           PPGL_PA_CS = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PA'|dataIn$Condition=='CS'),]},
           PPGL_PA_PHT = {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PA'|dataIn$Condition=='PHT'),]},
           PPGL_PA_CS_PHT =  {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PA'|dataIn$Condition=='CS'|dataIn$Condition=='PHT'),]},
           PPGL_PA_CS_PHT_HV =  {dataIn <- dataIn[which(dataIn$Condition=='PPGL'|dataIn$Condition=='PA'|dataIn$Condition=='CS'|dataIn$Condition=='PHT'|dataIn$Condition=='HV'),]},
           EHT_PHT_HV = { 
             dataIn$Condition <- as.character(dataIn$Condition)
             dataIn$Condition[dataIn$Condition =="PA" | dataIn$Condition =="PPGL" | dataIn$Condition =="CS"] <- "EHT"
             dataIn$Condition <- factor(dataIn$Condition)
             dataIn <- dataIn[which(dataIn$Condition=='EHT'|dataIn$Condition=='PHT'|dataIn$Condition=='HV'),]
           },
           EHT_PHT = {
             dataIn$Condition <- as.character(dataIn$Condition)
             dataIn$Condition[dataIn$Condition =="PA" | dataIn$Condition =="PPGL" | dataIn$Condition =="CS"] <- "EHT"
             dataIn$Condition <- factor(dataIn$Condition)
             dataIn <- dataIn[which(dataIn$Condition=='EHT'|dataIn$Condition=='PHT'),]
           },
           Gender = { 
             #do nothing
             },
           Centre = { 
             #do nothing
           },
           PAge = { 
             #do nothing
           },
           Intragroup = { 
             #do nothing
           },
           BAH_APA = { 
             dataIn <- dataIn[which(dataIn$Condition=='BAH'|dataIn$Condition=='APA'),]
             #do nothing [intra PA classification]
           }
           )
    
    # feature selection
    
    #drop the factors after all the data manipulation
    dataIn <- droplevels(dataIn)
   
          
    if (substr(as.character(scenarioIn$classifierMethod),1,5) =="RWeka") { 
          print('this is a RWeka classifier')
      
          # setup Rweka options
          if (regexpr("ASC", as.character(scenarioIn$classifierMethod))!=-1) {
            Weka_model <- make_Weka_classifier('weka.classifiers.meta.AttributeSelectedClassifier')
            
            if (regexpr("CFS", as.character(scenarioIn$classifierMethod))!=-1) {
              #E_value <- c('weka.attributeSelection.CfsSubsetEval') 
              E_value <- c('weka.attributeSelection.CfsSubsetEval -L -P 1 -E 1')
            }
            
            if (regexpr("Infogain", as.character(scenarioIn$classifierMethod))!=-1) {
              E_value <- c('weka.attributeSelection.InfoGainAttributeEval')
            }
            
            if (regexpr("_Boruta", as.character(scenarioIn$classifierMethod))!=-1) {
              # set Weka options to use all selected features i.e. Infogain + RankerAll
              # use Boruta later (just after data split and before Weka)
              E_value <- c('weka.attributeSelection.InfoGainAttributeEval')
              S_value <- c('weka.attributeSelection.Ranker')
            }
            
            
            if (regexpr("_Best", as.character(scenarioIn$classifierMethod))!=-1) {
              
              S_value <- c('weka.attributeSelection.GreedyStepwise')
              if (regexpr("_Best@F_", as.character(scenarioIn$classifierMethod))!=-1) {
              
                S_value <- paste(S_value,'-T -1.7976931348623157 -N -1 -num-slots 1', sep = ' ')
                  } else if (regexpr("_Best@B_", as.character(scenarioIn$classifierMethod))!=-1) {
                    S_value <- paste(S_value,'-D 0', sep = ' ')
                  } else if (regexpr("_Best@BD_", as.character(scenarioIn$classifierMethod))!=-1) {
                    S_value <- paste(S_value,'-D 2', sep = ' ')
                  }
              }
            
            if (regexpr("_RankerAll", as.character(scenarioIn$classifierMethod))!=-1) {
              S_value <- c('weka.attributeSelection.Ranker')
            }
            
            if (regexpr("_RF", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.trees.RandomForest')
            } else if (regexpr("_J48", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.trees.J48')
            } else if (regexpr("_LMT", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.trees.LMT')
            } else if (regexpr("_SL", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.functions.SimpleLogistic')
            } else if (regexpr("_SMO", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.functions.SMO')
            } else if (regexpr("_NB", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.bayes.NaiveBayes')
            } else if (regexpr("_IBk", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.lazy.IBk')
            } else if (regexpr("_Bag", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.meta.Bagging') # default with RepTree
            } else if (regexpr("_BN", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.bayes.BayesNet')
            } else if (regexpr("_LB", as.character(scenarioIn$classifierMethod))!=-1) {
              W_value <- c('weka.classifiers.meta.LogitBoost') # default with DS
            }
           
            optns <- Weka_control(E = list(E_value), S = list(S_value),
                                  W = list(W_value))
             
          }
      
          #set seed
          SEED <- scenarioIn$randomSeed
          set.seed(SEED,"L'Ecuyer") #set.seed(SEED,"L'Ecuyer") 
          
          
          switch(as.character(scenarioIn$validationType), 
                 DS_80_20 = {
                   # create a list of 80% of the rows in the original dataset we can use for training
                   if (scenarioIn$gensyn_flag == 1) {
                    
                    syn<-dataIn[dataIn$datatype==2,]
                    syn<-syn[,-which(colnames(syn) == 'datatype')]
                    orig<-dataIn[dataIn$datatype==1,]
                    orig<-orig[,-which(colnames(orig) == 'datatype')]
                    
                    spl = sample.split(orig$Condition , SplitRatio = 0.8)
                    dataTrain = subset(orig, spl==TRUE)
                    dataTest = subset(orig, spl==FALSE)
                    
                    dataTrain <- rbind(dataTrain,syn)
                    
                   }
                   else {
                   spl = sample.split(dataIn$Condition , SplitRatio = 0.8)
                   dataTrain = subset(dataIn, spl==TRUE)
                   dataTest = subset(dataIn, spl==FALSE) 
                   
                   dataTest$Condition <- droplevels(dataTest$Condition)
                   
                   if (length(levels(dataTest$Condition)) != length(levels(dataTrain$Condition))){
                     ## code added to check all classes in train are present in test else drop them in train
                   dataTrain <- dataTrain[which(dataTrain$Condition %in% levels(dataTest$Condition)),]
                   dataTrain$Condition <- droplevels(dataTrain$Condition)
                   warning_message <- paste0('All classes in Training set where not present in Test set, hence they were deleted from Training. Scenario number = ', scenarioIn$scenarioNumber, ", Random seed = ",
                                          scenarioIn$randomSeed, ", Condition = ", scenarioIn$conditionToClassify)
                   print(warning_message)
                   write(warning_message,file=paste0(scenarioIn$outDir,"Trained_model_warning.txt"),append=TRUE)
                   }
                  
                   }
                   
                   #if feature selection is Boruta
                   if (regexpr("_Boruta", as.character(scenarioIn$classifierMethod))!=-1) {
                     # set Weka options to use all selected features i.e. Infogain + RankerAll
                     # use Boruta later (just after data split and before RWeka)
                     
                     dir.create(file.path(scenarioIn$outDir, 'Boruta_files'), showWarnings = FALSE)
              
                     
                     temp_filename <- paste0(scenarioIn$outDir,"Boruta_files/Boruta_gen_data_", scenarioIn$conditionToClassify, "_", scenarioIn$randomSeed , "_output.RData")
                     if (file.exists(temp_filename)) {
                       #load data from saved pre feature selected files
                       load(temp_filename)
                     } else {
                       #process boruta if no file is found
                       boruta.train <- Boruta(Condition~., data = dataTrain, doTrace = 0)
                       final.boruta <- TentativeRoughFix(boruta.train)
                       
                       if (identical(getSelectedAttributes(final.boruta), character(0))) { # new addition for fixing scenarios which do not run)
                         fail_comment <- paste0('Boruta could not find appropriate features for Scenario number = ', scenarioIn$scenarioNumber, ", Random seed = ",
                                                scenarioIn$randomSeed, ", Condition = ", scenarioIn$conditionToClassify)
                         print(fail_comment)
                         write(fail_comment,file=paste0(scenarioIn$outDir,"Boruta_failed.txt"),append=TRUE)
                         
              
                          
                         out_data <- list('Accuracy'= NA, 'CorrectlyClassified' = NA,'IncorrectlyClassified' = NA,
                                          'TotalInstances' = NA, 'Balanced_Accuracy' = NA, 'Kappa' = NA, 'F1' = NA, 
                                          'Sensitivity' = NA, 'Specificity' = NA, 'Precision' = NA,'Recall' = NA, 
                                          'ComputeTime' = NA,'AUCScore' = NA) 
                         return(out_data)
                         
                       } else {
                       
                       dataTrain <- dataTrain[, names(dataTrain)[(names(dataTrain) %in% c(getSelectedAttributes(final.boruta),'Condition'))]]
                       dataTest <- dataTest[, names(dataTest)[(names(dataTest) %in% c(getSelectedAttributes(final.boruta),'Condition'))]]
                       save(dataTrain, dataTest, file= paste0(scenarioIn$outDir,"Boruta_files/Boruta_gen_data_", scenarioIn$conditionToClassify, "_", scenarioIn$randomSeed , "_output.RData"))
                       }
                       
                       
                       
                       
                     }
                     
                   }
                   
                   #start time tic
                   startTime <- Sys.time()
                   
                   trained_Weka_model <- Weka_model(dataTrain$Condition ~ ., dataTrain,control = optns)
                    
                   #end time toc
                   endTime <- Sys.time()
                   
                   #calculate elapsed time
                   elapsedTime <- endTime - startTime
                   
                   dataTest.pred <- predict(trained_Weka_model, newdata = dataTest, type = 'class') # add type = 'probability' for probabilities instead of labels 
                   dataTest.pred_prob <- predict(trained_Weka_model, newdata = dataTest, type = 'probability')
                   # can match the prediction with actual 
                   # table(dataTest$Condition, dataTest.pred)
                   # now add the metrics code from http://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
                   #  can also use the ML metrics package (esp. for MAE)
                   output <- confusionMatrix(dataTest.pred,dataTest$Condition, mode='everything')
              
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
                   
                   #calculate AUC
                   AUC <- auc(multiclass.roc(dataTest$Condition,dataTest.pred_prob))[1]
                   
                   if (is.na(Precision)) warning(paste0('No TP or FP in iteration: ',scenarioIn$scenarioNumber, ', hence Precision set to 0'))
                   #https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
                   Precision[is.na(Precision)] <- 0 # extreme cases when none classified correctly
                   
                   if (is.na(Recall)) warning(paste0('No TP or FN in iteration: ',scenarioIn$scenarioNumber, ', hence Recall set to 0'))
                   #https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
                   Recall[is.na(Recall)] <- 0 # extreme cases when none classified correctly
                   
                   if (is.na(F1)) warning(paste0('No TP or FN in iteration: ',scenarioIn$scenarioNumber, ', hence F1 set to 0'))
                   F1[is.na(F1)] <- 0 # extreme cases when none classified correctly
                   
              
                   
                   # buffer the java object before saving (else will be null)
                   .jcache(trained_Weka_model$classifier)
                   
                   # delete the bulky file before saving
                   rm(ensat_ht_data) 
                   
                   save(list = ls(all.names = TRUE), file= paste(scenarioIn$outDir,"Individual_Result_ScenarioNumber_", scenarioIn$scenarioNumber, "_output.RData", sep = '' ))
                   
                   out_data <- list('Accuracy'= Accuracy, 'CorrectlyClassified' = CorrectlyClassified,'IncorrectlyClassified' = IncorrectlyClassified,
                                    'TotalInstances' = TotalInstances, 'Balanced_Accuracy' = Balanced_Accuracy, 'Kappa' = Kappa, 'F1' = F1, 
                                    'Sensitivity' = Sensitivity, 'Specificity' = Specificity, 'Precision' = Precision,'Recall' = Recall, 
                                    'ComputeTime' = as.numeric(elapsedTime,units = "secs"),'AUCScore' = AUC)  
                   gc()
                   J("java.lang.Runtime")$getRuntime()$gc()
                   },
                 F_5 = {
                   # Run algorithms using 5-fold cross validation
                   #start time tic
                   startTime <- Sys.time()
                   
                   trained_Weka_model <- Weka_model(dataIn$Condition ~ ., dataIn,control = optns)
                   #end time toc
                   endTime <- Sys.time()
                   
                   #calculate elapsed time
                   elapsedTime <- endTime - startTime
                   
                   output <- evaluate_Weka_classifier(trained_Weka_model, dataIn, numFolds = 5, complexity = FALSE, seed = SEED, class = TRUE, type = c("probability")) 
              
                   
                   Accuracy <- output$details['pctCorrect']/100
                   CorrectlyClassified <- sum(diag(output$confusionMatrix)) 
                   IncorrectlyClassified <- sum(output$confusionMatrix) - CorrectlyClassified
                   TotalInstances <- sum(output$confusionMatrix)
                   WAvgPrecision <-sum(output$detailsClass[,'precision']*rowSums(output$confusionMatrix))/sum(rowSums(output$confusionMatrix))
                   WAvgRecall <- sum(output$detailsClass[,'recall']*rowSums(output$confusionMatrix))/sum(rowSums(output$confusionMatrix))
                   Kappa <- output$details['kappa']
                   WAvgFMeasure <- sum(output$detailsClass[,'fMeasure']*rowSums(output$confusionMatrix))/sum(rowSums(output$confusionMatrix))
                   
                   # buffer the java object before saving (else will be null)
                   .jcache(trained_Weka_model$classifier)
                   
                   # delete the bulky file before saving
                   rm(ensat_ht_data) 
                   
                   save(list = ls(all.names = TRUE), file= paste(scenarioIn$outDir,"Individual_Result_ScenarioNumber_", scenarioIn$scenarioNumber, "_output.RData",sep = ''))
                   
                   out_data <- list('Accuracy'= as.numeric(Accuracy), 'CorrectlyClassified' = CorrectlyClassified,'IncorrectlyClassified' = IncorrectlyClassified,'TotalInstances' = TotalInstances, 'WAvgPrecision' = WAvgPrecision,'WAvgRecall' = WAvgRecall,'WAvgFMeasure'= WAvgFMeasure,'Kappa'= Kappa,'ComputeTime' = as.numeric(elapsedTime,units = "secs"))  
                   gc()
                   J("java.lang.Runtime")$getRuntime()$gc()
                   
                   })
         
    }
    
    return(out_data)
}



