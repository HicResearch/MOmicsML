#############################################
# File:   scenarioResultAnalysis_Summary.r
#
# 
# Arguments: output folder path
#
# Outputs: performance metrics plots
#
# Usage: This script geneartes the various perfroamnce metrics plots for scenario rus.
#
#	Author: PS Reel

scenarioResultAnalysis_Summary <- function(outDirfinal) {
  library(ggplot2)
  library(plyr)
  library(stringr)
  
  scenarioData <- read.csv(paste(outDirfinal,'allScenariosSummaryResult.csv', sep=""),stringsAsFactors = TRUE)
  
  temp_frame <- data.frame(paste(rep("test_",5),c(1:5),sep=""))
  colnames(temp_frame)[1] <- "test"
  
  for(f in 1:length(levels(scenarioData$omicsData))){
    for(g in 1:length(levels(scenarioData$conditionToClassify))){
      
      scenarioData_subset <- scenarioData[which((scenarioData$omicsData %in% levels(scenarioData$omicsData)[f]) & (scenarioData$clinicalPhenotype %in% levels(scenarioData$clinicalPhenotype)) & (scenarioData$conditionToClassify %in% levels(scenarioData$conditionToClassify)[g]) & (scenarioData$validationType %in% 'DS_80_20')),]
      
      #check if the subset has valid results (i.e. no NA/NANs)
      if(sum(is.na(scenarioData_subset[,'Accuracy'])) != nrow(scenarioData_subset))
      {
      
      plot_out <- c("Accuracy", "Balanced_Accuracy","Kappa", "F1","Sensitivity","Specificity","Precision", "Recall","AUCScore")

      for (i in 1:length(plot_out)){
        print(i)
        a <- with(scenarioData_subset, reorder(scenarioData_subset$classifierMethod, scenarioData_subset[[plot_out[i]]], FUN = median))
        means_value <- aggregate(scenarioData_subset[[plot_out[i]]] ~  scenarioData_subset$classifierMethod, scenarioData_subset, mean)
        colnames(means_value) <-c("Algorithm",plot_out[i])

        if (grepl('Best', levels(scenarioData_subset$classifierMethod)[1])){
          label_data <- c("RWeka_ASC_CFS_Best@F_Bag" = "Bag", "RWeka_ASC_CFS_Best@F_BN" = "BN", "RWeka_ASC_CFS_Best@F_IBk" = "IBk", "RWeka_ASC_CFS_Best@F_J48" = "J48", "RWeka_ASC_CFS_Best@F_LB" = "LB",
                          "RWeka_ASC_CFS_Best@F_LMT" = "LMT", "RWeka_ASC_CFS_Best@F_NB" = "NB", "RWeka_ASC_CFS_Best@F_RF" = "RF",  "RWeka_ASC_CFS_Best@F_SL" ="SL",
                          "RWeka_ASC_CFS_Best@F_SMO" ="SMO")
        } else if (grepl('Boruta', levels(scenarioData_subset$classifierMethod)[1])) {
          label_data <- c("RWeka_ASC_Boruta_Bag" = "Bag", "RWeka_ASC_Boruta_BN" = "BN", "RWeka_ASC_Boruta_IBk" = "IBk", "RWeka_ASC_Boruta_J48" = "J48", "RWeka_ASC_Boruta_LB" = "LB",
                          "RWeka_ASC_Boruta_LMT" = "LMT", "RWeka_ASC_Boruta_NB" = "NB", "RWeka_ASC_Boruta_RF" = "RF",  "RWeka_ASC_Boruta_SL" ="SL",
                          "RWeka_ASC_Boruta_SMO" ="SMO")
        } else if (grepl('RankerAll', levels(scenarioData_subset$classifierMethod)[1])) {
          label_data <- c("RWeka_ASC_Infogain_RankerAll_Bag" = "Bag", "RWeka_ASC_Infogain_RankerAll_BN" = "BN", "RWeka_ASC_Infogain_RankerAll_IBk" = "IBk", "RWeka_ASC_Infogain_RankerAll_J48" = "J48", "RWeka_ASC_Infogain_RankerAll_LB" = "LB",
                          "RWeka_ASC_Infogain_RankerAll_LMT" = "LMT", "RWeka_ASC_Infogain_RankerAll_NB" = "NB", "RWeka_ASC_Infogain_RankerAll_RF" = "RF",  "RWeka_ASC_Infogain_RankerAll_SL" ="SL",
                          "RWeka_ASC_Infogain_RankerAll_SMO" ="SMO")
        }

        ggplot(scenarioData_subset, aes(a, scenarioData_subset[[plot_out[i]]])) +
          stat_boxplot(geom ='errorbar', colour = "#3366FF") +
          geom_boxplot(fill = "white", colour = "#3366FF", outlier.size = 0.8) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.margin = ggplot2::margin(10, 10, 10, 50)) +
          ggtitle(paste0(plot_out[i]," (RR - ", as.character(max(scenarioData_subset$randomSeed)),") \n ", levels(scenarioData$omicsData)[f], " for ", levels(scenarioData$conditionToClassify)[g])) +
          theme(plot.title = element_text(hjust=0.5,lineheight=.8, face="bold"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(size=12)) +
          coord_cartesian(ylim = c(0, 1)) +
          geom_hline(yintercept=(1/(str_count(levels(scenarioData$conditionToClassify)[g], "_")+1)), linetype="dashed", color = "red")+
          scale_x_discrete(labels=label_data)

        ggsave(paste(outDirfinal,levels(scenarioData$omicsData)[f], "-none-", levels(scenarioData$conditionToClassify)[g], "-DS_80_20_",plot_out[i],".pdf",sep = ""), width=4, height=4)
        ggsave(paste(outDirfinal,levels(scenarioData$omicsData)[f], "-none-", levels(scenarioData$conditionToClassify)[g], "-DS_80_20_",plot_out[i],".png",sep = ""), width=4, height=4)
      }
      
      }
      

      # find feature frequencies
      classifiers_used <- levels(scenarioData_subset$classifierMethod)
      
      for (i in 1:length(classifiers_used)){

        scenarioData_subset_classifier <- scenarioData[which((scenarioData$classifierMethod %in% classifiers_used[i]) & (scenarioData$omicsData %in% levels(scenarioData$omicsData)[f]) & (scenarioData$conditionToClassify %in% levels(scenarioData$conditionToClassify)[g])),]
        features <- scenarioData_subset_classifier[,which(colnames(scenarioData_subset_classifier)=="Feature_1"):ncol(scenarioData_subset_classifier)]
        features[] <- lapply(features, gsub, pattern='O1_', replacement='O1_')
        features[] <- lapply(features, gsub, pattern='O2_', replacement='O2_')
        features[] <- lapply(features, gsub, pattern='O3_', replacement='O3_')
        features[] <- lapply(features, gsub, pattern='O4_', replacement='O4_')
        features[] <- lapply(features, gsub, pattern='O5_', replacement='O5_')
        features[] <- lapply(features, gsub, pattern='O6_', replacement='O6_')
        features[] <- lapply(features, gsub, pattern='O7_', replacement='O7_')
        features <- stack (features)
        freqs <- ddply(features, .(values), summarise, count=length(values))
        freqs <- freqs[freqs$values!='', ]
        freqs$values <- factor(freqs$values,
                               levels=with(freqs,
                                           values[order(count, values, decreasing = TRUE)]))
        
        ggplot(freqs, aes(x=values,y=count)) +geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
          ggtitle(paste0(classifiers_used[i]," (RR - ", as.character(max(scenarioData_subset$randomSeed)),") ",levels(scenarioData$omicsData)[f]," ",levels(scenarioData$conditionToClassify)[g] ))+
          theme(plot.title = element_text(hjust=0.5,lineheight=.8, face="bold"),axis.text.x = element_text(size=8))
        ggsave(paste(outDirfinal,levels(scenarioData$omicsData)[f],"-none-",levels(scenarioData$conditionToClassify)[g],"-DS_80_20_",classifiers_used[i],".pdf",sep = ""), width=12, height=4)

        
        if (i==1){
          cutoff=0
          final_features <- data.frame(freqs[which(freqs$count >= (max(scenarioData_subset$randomSeed)*(cutoff/100)) ), 1])
          
          features_to_plot <- data.frame(freqs[which(freqs$count >= (max(scenarioData_subset$randomSeed)*(cutoff/100))),])
          features_to_plot <- features_to_plot[order(features_to_plot$count, decreasing = TRUE),]
          if (length(features_to_plot$values)!=0){
            pdf(paste(outDirfinal,levels(scenarioData$omicsData)[f],"-none-",levels(scenarioData$conditionToClassify)[g],"-DS_80_20_",classifiers_used[i],"_top50.pdf",sep = ""), width=12, height=4)
            a<-barplot(features_to_plot[,2])
            text(a[,1], -0.5, srt = 40, adj= 1, xpd = TRUE, labels = features_to_plot[,1], cex=1.2)
            dev.off()
          }
          
          colnames(features_to_plot)[1] <- paste(levels(scenarioData$conditionToClassify)[g],"F",sep = "_")
          colnames(features_to_plot)[2] <- paste(levels(scenarioData$conditionToClassify)[g],"C",sep = "_")
          temp_frame <-qpcR:::cbind.na(temp_frame,features_to_plot)
        }
        
      }
      
    }
  }
  temp_frame <- temp_frame[,-1]
  write.csv(temp_frame, file = paste(outDirfinal,"data_for_Venn_diagram_with_cutoff_", cutoff,".csv", sep=''), row.names = FALSE, na = '')
  print('Plots are successfully generated and saved')
  
}
