



BootstrapAnalysis <- function(Y, M, file.out="Bootstrap.out")
{

  year <- Y$YEAR[1]
  area <- Y$AREA[1]

  nb.data <- dim(Y)[1]
  mem1.estimates <- matrix(NA, ncol=4, nrow=M)
  sw <- rep(NA, M)
  sem1.estimates <- matrix(NA, ncol=4, nrow=M)
  mem2.estimates <- matrix(NA, ncol=4, nrow=M)
  sem2.estimates <- matrix(NA, ncol=4, nrow=M)

  lost <- 0
  for(i in 1:M)
  {
    ind.i<-sample(1:nb.data, nb.data, replace=TRUE)
    Ybis <- Y[ind.i,]
    while( sum(Ybis$N1+Ybis$N2)== 0 & lost < 1000)
    {
      ind.i<-sample(1:nb.data, nb.data, replace=TRUE)
      lost <- lost+1
    }
    if(lost==1000)
      {
        break("Not enough non zero data to perform extensive Bootstrap study")
      }
    mem1.estimates[i,] <- MEM.MLE(Ybis)
    sw[i] <- CPUE(Ybis)
    sem1.estimates[i,] <- SEM.MLE(Ybis)
    mem2.estimates[i,] <- MEM.MLE(Ybis, MEM=2)
    sem2.estimates[i,] <- SEM.MLE(Ybis)
  }
  
  ic.mem1=sort(mem1.estimates[,1])[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.sw=sort(sw)[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.sem1=sort(sem1.estimates)[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.mem2=sort(mem2.estimates[,1])[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.sem2=sort(sem2.estimates)[c(floor(M * 0.025), ceiling(M * 0.975))]

  estimates=matrix(c(mem1.estimates, mem2.estimates, sem1.estimates,  sem2.estimates, sw), ncol=17)
  estimates <- as.data.frame(estimates)
  names(estimates) <- c("mem1.lambda1", "mem1.lambda2", "mem1.p1", "mem1.p2",
                        "mem2.lambda1", "mem2.lambda2", "mem2.p1", "mem2.p2",
                        "sem1.lambda1", "sem1.lambda2", "sem1.lambdae", "sigma",
                        "sem2.lambda1", "sem2.lambda2", "sem2.lambdae", "sigma","sw")
  write.table(row.names=F, file=file.out, x=estimates, col.names=T)
  return(list(pt.estimates=estimates, IC.mem1=ic.mem1, IC.mem2=ic.mem2, IC.sem1=ic.sem1, IC.sem2=ic.sem2, IC.sw=ic.sw,  lost=lost) )
}

