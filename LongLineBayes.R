##############################################################
## Method 3  : Bayesian estimation 
##############################################################
## script to perform the Bayesian estimation of the complete exponetial model
## using BRugs package
##############################################################################
require(BRugs)


AnalysisBayes <- function(subdata, init=c(0.0001, 0.001), thin=50, burnin=2000, length.MC.chain=5000)
{
  data4Bugs <- list(x=matrix(c(subdata$N0, subdata$N1, subdata$N2), ncol=3), N= subdata$N, P=subdata$soaktime, L=nrow(subdata) )
  bugsData(data4Bugs, "data4Bugs.txt")
  bugsData(list(lambda=init), "init4Bugs.txt")
  modelCheck("E:/MarieEtienne/Longlines/OpenBugsCode/SeveralLonglinesExponentialModel.txt")
  modelData("data4Bugs.txt")
  modelCompile(numChains=1)
  modelInits("init4Bugs.txt")
  modelUpdate(burnin, thin=thin)
  samplesSet(c("lambda"))
  modelUpdate(length.MC.chain, thin=50)
  ok <- samplesStats("lambda")
  print(ok)
  return(list(
      lambda1=samplesSample("lambda[1]"), 
      lambda2=samplesSample("lambda[2]"),
	stat=ok
             )
        )
}

