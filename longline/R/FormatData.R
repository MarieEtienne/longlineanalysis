
# FormatDataForJags
#' FormatDataForJags prepare the data to be used with MEM model in jags
#' @param llData an object of class longline
#' @param parPrior a 2 components vector on the probability of espace pe~dbeta(parPrior[1], parPrior[2])
#' @param tau.lambda the precision used in the lognormal prior for all log.mu, log.lambda*F* parameters
#' @export
#' @return a list to be used as data.list in jags.model with model MEM
#

FormatDataForJags <- function(llData, parPrior=c(1,1), tau.lambda=0.5)
{
    Catch <- matrix(NA, nrow=llData$NData, ncol=4)
    Catch[,1]= llData$Nb
    Catch[,2]= llData$NT
    Catch[,3]= llData$NNT
    Catch[,4]= llData$Ne

  Fact1 <- as.numeric(as.factor(llData$Fact1))
  data.list=list(NFact1=llData$NFact1, NData=llData$NData,N=llData$N,
                 Catch=Catch, S=llData$S, Fact1=Fact1, a.moy=parPrior[1], b.moy=parPrior[1], 
                 log.lambda1F1=c(0, rep(NA, llData$NFact1-1)), log.lambda2F1=c(0, rep(NA, llData$NFact1-1)), tau.lambda=tau.lambda
  )
  
  
  if(with(llData, exists("Fact2")))
    data.list <- c(data.list, with(llData, list(Fact2=Fact2, NFact2=NFact2,log.lambda1F2=c(0, rep(NA, llData$NFact2-1)), log.lambda2F2=c(0, rep(NA, llData$NFact2-1))) ))
    
  return(data.list)
}

# FormatInitForJags
#' FormatInitForJags prepare the initial values to be used with MEM model in jags
#' @param llData an object of class longline
#' @param peConst should the probability of escape be considered as constant over Fact1, default is FALSE
#' @param  MEM stands for the version of MEM 1, if escape only arises from non target species, and 2 if all species have the same probability of escape
#' @export 
#' @return a list to be used as data.list in jags.model with model MEM
#

FormatInitForJags <- function(llData, peConst=F, MEM=1)
{
  estimate <- MEM.MLE(llData = llData, MEM=MEM)  
  ## using very small as initial values if no Target species observed 
  ## avoid issues with log(0) and no importance as it is just used for initialisation
  estimate$hat.lambda1[which(estimate$hat.lambda1<1e-13)] <- 1e-13
  NFact1 <- llData$NFact1
  NFact2 <- 1
  gg <- Reduce("rbind", strsplit(estimate$Group, "-"))
  if(nrow(estimate)==1){
     estimate$Fact1 <-gg[1]
  }else
  {
    estimate$Fact1 <-gg[,1]
    if(dim(gg)[2]>1 ){
      estimate$Fact2 <- gg[,2]
      NFact2 <- llData$NFact2    
    } 
  }
    
  

  if(NFact1==1 & NFact2<=1){
     formRight <- "1"
  }   else if (NFact1==1 & NFact2 > 1){
    formRight <- "as.factor(Fact2)"
  }   else if (NFact1>1 & NFact2 <= 1){
    formRight <- "as.factor(Fact1)"
  }   else if (NFact1>1 & NFact2 >1){
    formRight <- "as.factor(Fact1)+as.factor(Fact2)"
  }
  
  cat(paste0("with formula  ", formRight, "\n"))
  cat("\n *************************** \n")  
  
  
  coef.l1 <- coef(lm(formula=paste0("log(hat.lambda1)~", formRight), data=estimate))
  coef.l2 <- coef(lm(formula=paste0("log(hat.lambda2)~", formRight), data=estimate))
  log.mu1 <- coef.l1[1]
  log.mu2 <- coef.l2[1]
  log.lambda1F1 <- NA
  log.lambda2F1 <- NA
  if(NFact1>1){
    log.lambda1F1 <- c(log.lambda1F1, coef.l1[2:NFact1])
    log.lambda2F1 <- c(log.lambda2F1, coef.l2[2:NFact1])
  }
  
  log.lambda1F2 <- NULL
  log.lambda2F2 <- NULL
    if(NFact2>1){
      log.lambda1F2 <- c(NA, log.lambda1F2, coef.l1[(NFact1+1):(NFact1+NFact2-1)])
      log.lambda2F2 <- c(NA, log.lambda2F2, coef.l2[(NFact1+1):(NFact1+NFact2-1)])
    
  }
     pe <- mean(estimate$hat.pNT)
     if(!peConst)
       pe <- rep(pe, NFact1*NFact2)
     init.list <- list(log.lambda1F2=log.lambda1F2, log.lambda1F1=log.lambda1F1, log.lambda2F1=log.lambda2F1,log.lambda2F2=log.lambda2F2, pe=pe, log.mu1=log.mu1, log.mu2=log.mu2)
  init.list <- init.list[ ! sapply(init.list, is.null) ]
  
  return(init.list)
  
}
    