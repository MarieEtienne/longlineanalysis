# RunBayesEstimation
#' RunBayesEstimation runs the estimation of lambdaT and lambdaNT and produce abundance indices 
#' @param llData an object of class longline
#' @param MEM the version of the model, default is 1
#' @param nIter number of MCMC iterations kept, default is 2000
#' @param burnin  number of MCMC iterations dropped for the bruning period
#' @param filePath  the path where results and figure should be stored
#' @param parPrior a 2 components vector on the probability of espace pe~dbeta(parPrior[1], parPrior[2])
#' @param peConst should the probability of escape be considered as constant, default is FALSE prior for all log.mu*, log.lambda*F* parameters, default value=0.5
#' @param tau.lambda the precision used in the lognormal
#' @export
#' @return a list of abundance indice, and write some figures and tables in the directory specified in pathFile
#

RunBayesEstimation <- function(llData, MEM=1, nIter =2000, burnin = 500, filePath="", parPrior=c(1,1), peConst =F, tau.lambda=0.5)
{
  dataName <- as.list(sys.call() )$llData
  cat(paste0('Now running indices estimation on ', dataName, '\n')) 
  data.list  <-  FormatDataForJags(llData, tau.lambda=tau.lambda)
  init.list  <-  FormatInitForJags(llData, MEM = MEM, peConst = peConst)
  ## define Model name
  modelName <- paste0('MEM', MEM)
  if("Fact2" %in% names(llData)){
     modelName <- paste0(modelName, 'AdjustFact')
  }
  if(peConst)
    modelName <- paste0(modelName,'ConstantPe')
  
  modelName <- paste0(modelName,'.txt')
     
print(modelName)     
  jagsRun <- jags.model(file=system.file(package='longline', 'jags', modelName), data = data.list, inits = init.list, n.adapt = 100, n.chains = 1)
  update(jagsRun, burnin)
  out.samples <- jags.samples(jagsRun, c("log.mu1","log.mu2", "log.lambda1F1", "log.lambda2F1", "pe"), n.iter=nIter)

postIndices <- PostIndicesPlot(llData, out.samples, who="TF1", yLim = range(c(0)),xLab = levels(llData$Fact1), fileOut = file.path(filePath, paste0(dataName, '-MEM=',MEM, 'peConst=', peConst, '.pdf')) )

  ### save indices
  res <- Reduce('rbind',tapply(postIndices$index, postIndices$Fact1, 
                               function(d) {
                                 c(summary(d), sd(d), sd(d)/mean(d))
                               }))
  colnames(res) <- c(names(summary(rnorm(10))), "sd", "cv")
  rownames(res) <- rep("", nrow(res))
  res <- cbind(data.frame(Fact1=levels(llData$Fact1)), res)
  
  write.table(file = file.path(filePath,  paste0(dataName,'-MEM=', MEM,'-peConst=', peConst, '.csv')), res, sep=";", row.names=F)
  save(file = file.path(filePath,  paste0(dataName,'-MEM=', MEM,'-peConst=', peConst, '.Rd')), list="out.samples" )
  return(postIndices)
}

# postMean <- function(longline, outmcmc,  peConst=F)
# {
#   NYear <- longline$NYear
#   NArea <- longline$NArea
#   outl1 <- apply(outmcmc$lambda1,1, mean)
#   outl2 <- apply(outmcmc$lambda2,1, mean)
#   outpe <- apply(outmcmc$pe, 1, mean)
#   estimate <- data.frame(Year=rep(NA, NYear*NArea), Area=rep(NA, NYear*NArea),
#                          lambdaT=rep(NA, NYear*NArea), lambdaNT=rep(NA, NYear*NArea), 
#                          pe=rep(NA, NYear*NArea))
#   
#   yearunique = as.numeric(levels(as.factor(longline$Year)))
#   areaunique = as.numeric(levels(as.factor(longline$Area)))
#   compt=1
# 
#   
# #   for( a in 1:NArea){
# #     provlambdaT= data.frame(Year=NULL, Value=NULL)
# #     provlambdaNT= data.frame(Year=NULL, Value=NULL)
# #     for( y in 1:NYear)
# #     {
# #       if(sum(longline$Year==yearunique[y] & longline$Area==areaunique[a])>1)
# #       {
# #         n <- length(outmcmc$lambda1[1,,])
# #         p1 <- data.frame(Year=rep(yearunique[y],n), Value=outmcmc$lambda1[(a-1)*NYear+y,,])
# #         p2 <- data.frame(Year=rep(yearunique[y],n), Value=outmcmc$lambda2[(a-1)*NYear+y,,])
# #         provlambdaT = rbind(provlambdaT, p1)
# #         provlambdaNT = rbind(provlambdaNT, p2)
# #       }
# #     }
#     
#     
# for(a in 1:NArea)   
#   {
#   for(y in 1:NYear)
#     {
#       ind <- which(longline$Year==yearunique[y] & longline$Area==areaunique[a] )
#       if(length(ind>0)){
#         estimate[compt,1] <- yearunique[y]
#         estimate[compt,2] <- areaunique[a] 
#         estimate[compt,3] <- mean(outmcmc$lambda1[(a-1)*NYear + y ,,])
#         estimate[compt,4] <- mean(outmcmc$lambda2[(a-1)*NYear+y,,])
#         estimate[compt,5] <- ifelse( peConst,mean(outmcmc$pe) , mean(outmcmc$pe[(a-1)*NYear+y,,]))
#                                  #outl1[(a-1)*NYear+y], outl2[(a-1)*NYear+y], outpe[(a-1)*NYear+y])
#         compt=compt+1
#       }
#     }
#   }
#   if(compt<=nrow(estimate))
#   {
#     estimate <- estimate[-(compt:nrow(estimate)),]
#   }
#   return(estimate)
# 
# }
# 
# postSummary <- function(longline, outmcmc, peConst=F)
# {
#   NYear <- longline$NYear
#   NArea <- longline$NArea
#   Year  <- rep(NA, NYear*NArea)
#   Area  <- rep(NA, NYear*NArea)
#   outl1 <- apply(outmcmc$lambda1,1, summary)
#   outl2 <- apply(outmcmc$lambda2,1, summary)
#   outpe <- apply(outmcmc$pe, 1, summary)
#   nrows.tot <- NYear*NArea
#   estimate <- list(lambdaT=matrix(NA, nrow=NYear*NArea, ncol=8), lambdaNT=matrix(NA, nrow=NYear*NArea, ncol=8), 
#                          pe=matrix(NA, nrow=NYear*NArea, ncol=8)
#                    )
#   
#   colnames(estimate$lambdaT)=c(unlist(strsplit(outl1[,1],":"))[seq(1,11,2)], "sd", "cv")
#   colnames(estimate$lambdaNT)=c(unlist(strsplit(outl2[,1 ],":"))[seq(1,11,2)], "sd", "cv")
#   colnames(estimate$pe)=c(unlist(strsplit(outl1[,1],":"))[seq(1,11,2)],"sd", "cv")
#   yearunique = as.numeric(levels(as.factor(longline$Year)))
#   areaunique = as.numeric(levels(as.factor(longline$Area)))
#   compt=1
#   for(a in 1:NArea)
#     {
#     for(y in 1:NYear)
#     {
#       outl1num=as.numeric(unlist(strsplit(outl1[,(a-1)*NYear+y],":"))[seq(2,12,2)])
#       outl2num=as.numeric(unlist(strsplit(outl2[,(a-1)*NYear+y],":"))[seq(2,12,2)])
#       if(! peConst)
#         outpenum=as.numeric(unlist(strsplit(outpe[,(a-1)*NYear+y],":"))[seq(2,12,2)])
#       else
#         outpenum=as.numeric(unlist(strsplit(outpe,":"))[seq(2,12,2)])
#       
#       if(sum(longline$Year==yearunique[y] & longline$Area==areaunique[a])>1)
#       {
#         Year[compt] <- yearunique[y]
#         Area[compt] <- areaunique[a] 
#         estimate$lambdaT[compt,1:7] <- c(outl1num,sd(outmcmc$lambda1[(a-1)*NYear+y,,]))
#         estimate$lambdaNT[compt,1:7] <- c(outl2num, sd(outmcmc$lambda2[(a-1)*NYear+y,,]))
#         estimate$lambdaNT[compt,8] <- estimate$lambdaNT[compt,7] / estimate$lambdaNT[compt,4]
#         estimate$lambdaT[compt,8] <- estimate$lambdaT[compt,7] / estimate$lambdaT[compt,4]
#           if(!peConst){
#         estimate$pe[compt,1:7] <- c(outpenum, sd(outmcmc$pe[(a-1)*NYear+y,,]))
#           }else{
#             estimate$pe[compt,1:7] <- c(outpenum, sd(outmcmc$pe))
#           }
#         estimate$pe[compt,8] <- estimate$pe[compt,7] / estimate$pe[compt,4]
#               
#         compt=compt+1
#       }
#     }
#   }
#     
#     
#   if(compt<=nrows.tot)
#   {
#     Year <- Year[-(compt:nrows.tot )]
#     Area <- Area[-(compt:nrows.tot )]
#     estimate$lambdaT <- estimate$lambdaT[-(compt:nrows.tot ),]
#     estimate$lambdaNT <- estimate$lambdaNT[-(compt:nrows.tot ),]
#     estimate$pe <- estimate$pe[-(compt:nrows.tot ),]
#     
#   }
#   rownames(estimate$lambdaT)= (paste(Area, Year, sep="-"))
#   rownames(estimate$lambdaNT)=rownames(estimate$lambdaT)
#   return(estimate)
#   
# }