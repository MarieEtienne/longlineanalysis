# PostIndicesPlot
#' PostIndicesPlot is designed to produce plots for the posterior distribution of abundance indices for Target and non Target species, it may be used to produce posteriori boxplot of other quantitties
#' @param llData an object of class longline
#' @param outSamples the results of RunBayesEstimation on llData
#' @param who the quantity to be plotted if TF1 (resp NTF1 ) it produces the posterior distribution for target (resp Non Targe) species for every level of factor F1, Who may be any other posterior quantity of the model (log.mu1, log.mu2, log.lambda1F1, log.lambda2.F1, pe)
#' @param fileOut the file to save the plots
#' @param yLim optionnal the limits of the Y axis, default is null and the scales will be determined automatically
#' @param xLab optionnal the names of xlabel
#' @import scales ggplot2 
#' @export
#' @return no return value

PostIndicesPlot<- function(llData, outSamples, who="TF1", fileOut="plot.pdf", yLim=NULL, xLab=NULL, 
                           complete.missing.year=F){

  if(complete.missing.year){
    ypres <- unique(as.numeric(as.character(llData$Fact1)))
    miss <- setdiff(min(ypres  ):max(ypres), ypres)
  }
  if(who=="TF1"){
  out <- exp(sweep(matrix(outSamples$log.lambda1F1[,,], ncol=dim(outSamples$log.lambda1F1)[1], byrow=T), MARGIN = 1, FUN = "+", STATS = as.numeric(outSamples$log.mu1)))

  xLab = levels(llData$Fact1)
  } else {
    if ( who=="NTF1")
    {
      out <- exp(sweep(matrix(outSamples$log.lambda2F1[,,], ncol=dim(outSamples$log.lambda2F1)[1], byrow=T), MARGIN = 1, FUN = "+", STATS = as.numeric(outSamples$log.mu2)))
      xLab = levels(llData$Fact1)
    } else {
       if( who %in% names(outSamples)){
         out <- matrix(outSamples[[who]][,,], ncol=dim(outSamples[[who]])[1], byrow=T)
         xLab = levels(llData$Fact1)
       } else {
           cat(paste0(who, " is not a valid name for posterior plots. Possible values are ", names(outSamples)," \n"))
         }
    }
  }
  yLim <- range(c(yLim, range(out)))
  out2 <- data.frame(index=as.numeric(out), 
                               Fact1=as.factor(rep(1:dim(out)[2], each = dim(out)[1])))
  out2$Fact1 <- levels(llData$Fact1)[out2$Fact1]
    
  if(complete.missing.year){
    out2 <- rbind(out2, data.frame(Fact1=as.character(miss), index=rep(NA, length(miss))))
    out2$Fact1 <- as.factor(as.numeric(out2$Fact1))
    xLab = levels(out2$Fact1)
  }
  
  p <- ggplot(out2, aes(x=Fact1, y=index))
    p <- p + geom_boxplot() + scale_y_continuous(label=scientific, limits=yLim) + scale_x_discrete(labels=xLab  )
    p
  print(fileOut)
  ggsave(plot = p, filename = fileOut,width = 20, height =15, units = "cm", dpi = 300)
  if(complete.missing.year){
   out2 <- out2[!is.na(out2$index),]
   levels(out2$Fact1) <- droplevels(out2$Fact1)
     }
  return(out2)
}

 