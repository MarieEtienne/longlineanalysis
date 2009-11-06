######################################################
##   SimStudy launches a bunch of simu according
##   to the config file file.in
######################################################


SimStudy <- function( file.in )
  {
    par.simu <- ReadParameters(file.in)
    n.simu <- par.simu$n.simu
    
    file.out <- paste(strsplit(file.in,".init"), ".out", sep="")
    
    write.table(file=file.out, x=par.simu, sep=",", row.names=FALSE)
    Add2File(x=c("","mle.lambda1", "mle.lambda2", "bayes.lambda1", "bayes.lambda2", "nlr.lambda1", "nlr.lambda2","nlr.lambdae"), file=file.out)
    
    for( i in 1:n.simu)
      {
        Y <- DrawSampleswEscape(par.simu=par.simu)
        bayes <- AnalysisBayes(Y)
        mle <- AnalysisMLE(Y)
        nlr <- AnalysisNLR2(Y)
        sw <- computeSweptArea(Y)
        
        Add2File(x=c("Empty", mle$lambda, bayes$lambda, nlr$lambda, sw), file=file.out)  #Add to the results file
        
        ##empty hooks are considered with all other species
        Y$N <- Y$N
        Y$N2 <- Y$N2+Y$Ne
        Y$Ne <- rep(0,nrow(Y))
        
        nlr.empty <- AnalysisNLR2(Y)
        mle.empty <- AnalysisMLE(Y)
        bayes.empty <- AnalysisBayes(Y)
        
        Add2File(x=c("NoEmpty", mle.empty$lambda, bayes.empty$lambda, nlr.empty$lambda, sw), file=file.out)  #Add to the results file
        
        
      }
  }


#### from a real pattern
SimStudy2 <- function( lambda1, lambda2, Y, M, type="uniform", percent=0.00, pref=1, file.out="output.out") 
  {
    
    
    write.table(file=file.out, x=list(lambda1, lambda2, type, percent, pref) , sep=",", row.names=FALSE)
    Add2File(x=c("","mle.lambda1", "mle.lambda2", "nlr.lambda1", "nlr.lambda2","nlr.lambdae"), file=file.out)
    
    for( i in 1:M)
      {
        Y <- DrawSampleswEscape2(lambda1, lambda2, Y, type, percent, pref)
        mle <- AnalysisMLE(Y)
        nlr <- AnalysisNLR2(Y)
        cpue <- computeSweptArea(Y)
        
        Add2File(x=c("Empty", mle$lambda, nlr$lambda, sw), file=file.out)  #Add to the results file
        
        ##empty hooks are considered with all other species
        Y$N <- Y$N
        Y$N2 <- Y$N2+Y$Ne
        Y$Ne <- rep(0,nrow(Y))
        
        nlr.empty <- AnalysisNLR2(Y)
        mle.empty <- AnalysisMLE(Y)
        Add2File(x=c("NoEmpty", mle.empty$lambda, nlr.empty$lambda, sw), file=file.out)  #Add to the results file
        
        
      }
  }

