######################################################
##   SimStudy launches a bunch of simu according
##   to the config file file.in
##   forget.empty=0, empty hooks are taken into acoount in some way
##   forget.empty=1 empty hooks are just ignored, N is still N
##   forget.empty=2 empty hooks are removed and N=N-Ne
######################################################


SimStudy <- function( file.in, forget.empty=0 )
  {

    par.simu <- ReadParameters(file.in)
    n.simu <- par.simu$n.simu
    
    file.out <- paste(strsplit(file.in,".init"), ".out", sep="")
    
    write.table(file=file.out, x=par.simu, sep=",", row.names=FALSE)
    Add2File(x=c("","mle.lambda1", "mle.lambda2", "nlr.lambda1", "nlr.lambda2","nlr.lambdae", "sw"), file=file.out)
    
    for( i in 1:n.simu)
      {
      if(i%%50==0) print(i)
        Y <- DrawSampleswEscape(par.simu=par.simu)
        ## if empty hooks are forgotten
        if(forget.empty==1)
          {
            Y$Nb <- Y$Nb + Y$Ne
            Y$Ne <- 0 * Y$Ne
            
          }else if(forget.empty==2) 
        {
          Y$N <- Y$N- Y$Ne
          Y$Ne <- 0 * Y$Ne
        }                
#        bayes <- AnalysisBayes(Y)

         nlr <- AnalysisNLR(Y)
         mle <- AnalysisMLE(Y)
         sw <- computeSweptArea(Y)
        Add2File(x=c("Empty", mle$lambda, nlr$lambda, sw), file=file.out)  #Add to the results file
        
#       ##empty hooks are considered with all other species
#        Y$N <- Y$N
#        Y$N2 <- Y$N2+Y$Ne
#        Y$Ne <- rep(0,nrow(Y))
#        
#        nlr.empty <- AnalysisNLR(Y)
#        mle.empty <- AnalysisMLE(Y)
# #       bayes.empty <- AnalysisBayes(Y)
#        
#        Add2File(x=c("NoEmpty", mle.empty$lambda, nlr.empty$lambda, sw), file=file.out)  #Add to the results file
#        
        
      }
  }


#### from a real pattern
SimStudy2 <- function( lambda1, lambda2, Y, M, type="uniform", percent=0.00, pref=1, file.out="Simu.out") 
  {
    
  nb.data <- dim(Y)[1]
  mle.estimates <- matrix(NA, ncol=2, nrow=M)
  swept <- rep(NA, M)
  nlr.estimates <- matrix(NA, ncol=2, nrow=M)
  mle.estimates.no <- matrix(NA, ncol=2, nrow=M)
  nlr.estimates.no <- matrix(NA, ncol=2, nrow=M)
if(lambda1>0)
  {
    for( i in 1:M)
      {
        Y <- DrawSampleswEscape2(lambda1, lambda2, Y, type, percent, pref)
        mle.estimates[i,] <- AnalysisMLE(Y)$lambda
        nlr.estimates[i,] <- AnalysisNLR(Y)$lambda[c(1,2)]
        swept[i] <- computeSweptArea(Y)
        
        ##empty hooks are considered with all other species
        Ybis <- Y
        Ybis$N2 <- Ybis$N2+Ybis$Ne
        Ybis$Ne <- rep(0,nrow(Ybis))
        
        
        mle.estimates.no[i,] <- AnalysisMLE(Ybis)$lambda
        nlr.estimates.no[i,] <- AnalysisNLR(Ybis)$lambda[c(1,2)]
        
        
      }
  }
    estimates=matrix(c(mle.estimates, mle.estimates.no, nlr.estimates,  nlr.estimates.no, swept), ncol=9)
    estimates <- as.data.frame(estimates)
    names(estimates) <- c("MEM2.lambda1","MEM2.lambda2","MEM1.lambda1","MEM1.lambda2","SEM2.lambda1","SEM2.lambda2","SEM1.lambda1","SEM1.lambda2", "swept")
    write.table(row.names=F, file=file.out, x=estimates, col.names=T)

  }

