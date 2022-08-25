#############################################
# File:   saveConfusionMatrix.r
#
# 
# Arguments: results output folder path, scenario number and the confusion matrix 
#
# Outputs: confusion matrix in png format
#
# Usage: This script saves confusion matrix from the scenario runs in png image format.
#
#	Author: PS Reel

saveConfusionMatrix <- function(outDirfinal,i,x1) {
  library(ggplot2)
  if (i==1){
    dir.create(file.path(outDirfinal, 'ConfusionMatrix'))
  }
  
  png(paste0(outDirfinal,"/ConfusionMatrix/","CM_",i,".png"))
      
  x1 <- t(x1)
  opar <- par(mar=c(6.1, 7.1, 7, 6))
  x <- unclass(x1)
  
  
  diag(x) <- -diag(x)
  x <- -x
  image(1:ncol(x), 1:ncol(x), xlab = '', ylab='',
        t(apply(x, 2, rev)), xaxt='n', yaxt='n', 
        col=colorRampPalette(c('#e5f5f9','#2c7fb8'))(3),
        #col=colorRampPalette(c(hsv(h = 0, s = 0.9, v = 0.9, alpha = 1), 
        #                       hsv(h = 0, s = 0, v = 0.9, alpha = 1), 
        #                       hsv(h = 2/6, s = 0.9, v = 0.9, alpha = 1)))(16),
        zlim=c(min((x)), max((x))), cex=1.4) #intally set to -10,10
  axis(3, at=1:ncol(x), labels=colnames(x), las=2,cex.axis=0.8, cex.axis=1.4)
  axis(2, at=ncol(x):1, labels=colnames(x), las=1, cex.axis=0.8, cex.axis=1.4)
  axis(1, at=1:ncol(x), labels=colSums(x1), las=2,cex.axis=0.8, cex.axis=1.4)
  axis(4, at=ncol(x):1, labels=rowSums(x1), las=1, cex.axis=0.8, cex.axis=1.4)
  title(xlab='Actual Class', cex.lab=1.4)
  title(ylab='Predicted Class', line=6, cex.lab=1.4)
  abline(h = 0:ncol(x) + 0.5, col = 'gray')
  abline(v = 0:ncol(x) + 0.5, col = 'gray')
    
  text(rep(1:nrow(x1), each=nrow(x1)),nrow(x1):1, 
       labels = as.character(paste0(replicate(nrow(x1)*nrow(x1),'') ,c(x1))),cex = 1.4)
  box(lwd=2)
  
  
  title(main = paste0('Confusion Matrix', ' (', i,')'),line=6, cex.lab=1.4)
    
  dev.off()
  par(opar)
  
}