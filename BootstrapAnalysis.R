



BootstrapAnalysis <- function(Y, M, file.out="Bootstrap.out")
{

  year <- Y$YEAR[1]
  area <- Y$AREA[1]

  nb.data <- dim(Y)[1]
  mle.estimates <- matrix(NA, ncol=2, nrow=M)
  swept <- rep(NA, M)
  nlr.estimates <- matrix(NA, ncol=2, nrow=M)
  mle.estimates.no <- matrix(NA, ncol=2, nrow=M)
  nlr.estimates.no <- matrix(NA, ncol=2, nrow=M)

  lost <- 0
  for(i in 1:M)
  {
    ind.i<-sample(1:nb.data, nb.data, replace=TRUE)
    Ybis <- Y[ind.i,]
    while( sum(Ybis$N1)== 0 & lost < 1000)
    {
      ind.i<-sample(1:nb.data, nb.data, replace=TRUE)
      lost <- lost+1
    }
    if(lost==1000)
      {
        break("Not enough non zero data to perform extensive Bootstrap study")
      }
    mle.estimates[i,] <- AnalysisMLE(Ybis)$lambda
    swept[i] <- computeSweptArea(Ybis)
    nlr.estimates[i,] <- AnalysisNLR2(Ybis)$lambda[c(1,2)]

    Ybis$N2 <- Ybis$N2+Ybis$Ne
    Ybis$Ne <- 0*Ybis$Ne

    mle.estimates.no[i,] <- AnalysisMLE(Ybis)$lambda
    nlr.estimates.no[i,] <- AnalysisNLR2(Ybis)$lambda[c(1,2)]
    
  }
  ic=sort(mle.estimates[,1])[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.swept=sort(swept)[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.nlr=sort(nlr.estimates)[c(floor(M * 0.025), ceiling(M * 0.975))]
  estimates=matrix(c(mle.estimates, mle.estimates.no, nlr.estimates,  nlr.estimates.no, swept), ncol=9)
  estimates <- as.data.frame(estimates)
  names(estimates) <- c("MEM2.lambda1","MEM2.lambda2","MEM1.lambda1","MEM1.lambda2","SEM2.lambda1","SEM2.lambda2","SEM1.lambda1","SEM1.lambda2", "swept")
  write.table(row.names=F, file=file.out, x=estimates, col.names=T)
  return(list(pt.estimates=estimates, IC.mle=ic, IC.swept=ic.swept, IC.nlr=ic.nlr,  lost=lost) )
}

