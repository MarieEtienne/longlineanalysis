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

    Add2File(x=c("mem1.lambda1", "mem1.lambda2", "mem1.p1", "mem1.p2",
                 "mem1.lambda1.g", "mem1.lambda2.g", "mem1.p1.g", "mem1.p2.g",
                        "mem2.lambda1", "mem2.lambda2", "mem2.p1", "mem2.p2",
                        "mem2.lambda1.g", "mem2.lambda2.g", "mem2.p1.g", "mem2.p2.g",
                        "sem1.lambda1", "sem1.lambda2", "sem1.lambdae", "sigma",
                        "sem2.lambda1", "sem2.lambda2", "sem2.lambdae", "sigma","sw"), file=file.out)
    
    for( i in 1:n.simu)
      {
        Y <- DrawSampleswEscape(par.simu=par.simu)
        ## if empty hooks are forgotten
        if(forget.empty==1)
          {
            Y$Ne <- 0 * Y$Ne
          }else if(forget.empty==2) 
        {
          Y$N <- Y$N- Y$Ne
          Y$Ne <- 0 * Y$Ne
        }

        ##bayes <- AnalysisBayes(Y)
        mem1 <- MEM.MLE(Y)
        mem1.g <- MEM.MLE(Y, sameP=F)
        sem1 <- SEM.MLE(Y)
        mem2 <- MEM.MLE(Y, MEM=2)
        mem2.g <- MEM.MLE(Y, MEM=2, sameP=F)
        sem2 <- SEM.MLE(Y, SEM=2)
        sw   <- CPUE(Y) 
        
        Add2File(x=c(mem1, mem1.g, mem2.g, mem2, sem1, sem2, sw), file=file.out)  #Add to the results file        
      }
  }


#### from a real pattern
SimStudy2 <- function( lambda1, lambda2, Y, M, type="uniform", percent=0.00, pref=1, file.out="Simu.out") 
  {
    
  nb.data <- dim(Y)[1]
  mem1.estimates <- matrix(NA, ncol=4, nrow=M)
  mem2.estimates <- matrix(NA, ncol=4, nrow=M)
  sem1.estimates <- matrix(NA, ncol=4, nrow=M)
  sem2.estimates <- matrix(NA, ncol=4, nrow=M)
  sw <- rep(NA, M)
  if(lambda1>0)
  {
    for( i in 1:M)
      {
        Y <- DrawSampleswEscape2(lambda1, lambda2, Y, type, percent, pref)
        mem1.estimates[i,] <- MEM.MLE(Y)
        sem1.estimates[i,] <- SEM.MLE(Y)
        mem2.estimates[i,] <- MEM.MLE(Y, MEM=2)
        sem2.estimates[i,] <- SEM.MLE(Y, SEM=2)
        sw[i] <- CPUE(Y)
      }
  }
    estimates=matrix(c(mem1, mem2, sem1,sem2, sw), ncol=17)
    estimates <- as.data.frame(estimates)
    names(estimates) <- c("mem1.lambda1", "mem1.lambda2", "mem1.p1", "mem1.p2",
                        "mem2.lambda1", "mem2.lambda2", "mem2.p1", "mem2.p2",
                        "sem1.lambda1", "sem1.lambda2", "sem1.lambdae", "sigma",
                        "sem2.lambda1", "sem2.lambda2", "sem2.lambdae", "sigma","sw")
    write.table(row.names=F, file=file.out, x=estimates, col.names=T)

  }

