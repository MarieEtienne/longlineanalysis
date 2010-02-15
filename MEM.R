##############################################################
## Method 2  : ML estimation 
##############################################################
## script to perform the ML estimation of the complete exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## subdata the dataframe which stored the data to be analyzed
##############################################################################

AnalysisMLE <- function(subdata, tune1=3, sameP=TRUE, MEM=1) 
{
  ## initial values for minimization algorithm
  if(!is.CLonglineData(subdata))
  {
    stop("Argument 1 of function AnalysisMLE is not of class CLOngloneData")
  } else if( sum( subdata$N1+subdata$N2) == 0)
  {
    return(list( lambda = c( NA,NA, NA,NA ) ) ) ## no possible estimation
  } else if(sameP)     
  {
    if(MEM==1) 
    {
      return(Mem1MleSameP(subdata))
    }else if(MEM==2)
    {
      return(Mem2MleSameP(subdata))
    }             
  } else
  {
    if(MEM==1) 
    {
      return(Mem1MleGeneral(subdata))
    } else
    {
      return(Mem1MleGeneral(subdata))
    }
  }
}


################################################
## Estimators for MEM1 with shared soaktime P
###############################################
Mem1MleSameP <- function(subdata)
{
  with(subdata,
    {
      hat.lambda1=sum( N1 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.lambda2=sum( N2 + Ne) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.p2 = sum(Ne)/ sum(Ne+N2)
      hat.p1 =0
      return(c(hat.lambda1, hat.lambda2, hat.p1, hat.p2))
    }
  )
}

################################################
## Estimators for MEM2 with shared soaktime P
###############################################
Mem2MleSameP <- function(subdata)
{
  with(subdata,
    {
      hat.lambda1=sum( N1 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.lambda2=sum( N2 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.p2 = sum( Ne )/ sum( N - Nb )
      hat.p1 =  hat.p2 
      return(c(hat.lambda1, hat.lambda2, hat.p1, hat.p2))
    }
  )
}

################################################
## Estimators for MEM1 general case
###############################################
Mem1MleGeneral <- function(subdata)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( N2 + Ne ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + N2)
      theta.mem1 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p) ) )
      optim.mem1 = optim(theta.mem1, LogLike.Mem1, method="BFGS", control=list(fnscale=-1),  subdata=subdata)
      par.mem1   = optim.mem1$par 
      return( c(exp(par.mem1[1:2]), 0, exp(par.mem1[3]) / (1 + exp(par.mem1[3]) )) )
    }
  )
}


################################################
## Estimators for MEM2 general case
###############################################
Mem2MleGeneral <- function(subdata)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( N2 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + N2 + N1 )
      theta.mem2 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      
      optim.mem2 = optim(theta.mem2, LogLike.Mem2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      par.mem2   = optim.mem2$par 
      return( c(par.mem2, par.mem2[3]) )
    }
  )
}



######################################################
##   General LogLikelihood
######################################################
LogLike <- function(theta,  subdata)
{
  ## theta = c(log(lambda1), log(lambda2), logit(p1), logit(p2))
  lambda1 = exp(theta[1])
  lambda2 = exp(theta[2])
  p1      = exp(theta[3]) / ( 1 + exp(theta[3]) )
  p2      = exp(theta[4]) / ( 1 + exp(theta[4]) )
  
  with(subdata,
    {
      lambda=lambda1+lambda2
      return (
        sum( lgamma(N+1)-lgamma(Nb+1)-lgamma(N1+1)-lgamma(N2+1)-lgamma(Ne+1) ) - sum(lambda * P * Nb) +
        sum( (N-Nb) * log (1 - exp( -lambda * P) ) ) + sum(N1 * log(  lambda1 / lambda * ( 1-p1) ) )+ 
        sum( N2 * log(  lambda2 / lambda * ( 1-p2) ) ) + sum ( Ne * log(  (lambda1 *p1 + lambda2 * p2) / lambda ) )
        )
    }
  )
}

######################################################
##   MEM1 LogLikelihood
######################################################
LogLike.Mem1 <- function(theta.mem1,  subdata)
{
  ## theta.mem1 = c(log(lambda1), log(lambda2), logit(p2))
  return( LogLike( c(theta.mem1[1:2], -Inf, theta.mem1[3]), subdata ) )
}

######################################################
##   MEM2 LogLikelihood
######################################################
LogLike.Mem2 <- function(theta.mem2,  subdata)
{
  ## theta.mem1 = c(lambda1, lambda2, p)
  return( LogLike( c(theta.mem2[1:2], theta.mem2[3], theta.mem2[3] ), subdata ) )
}


######################################################
##   MEM1 AIC
######################################################
AIC.Mem1 <- function(subdata)
{
   with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( N2 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + N2 )
      theta.mem1 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      
      optim.mem1 = optim(theta.mem1, LogLike.Mem1, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem2$value + 6)
    }
    )
}

######################################################
##   MEM2 AIC
######################################################
AIC.Mem2 <- function(subdata)
{
   with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( N2 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + N2 + N1 )
      theta.mem2 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      
      optim.mem2 = optim(theta.mem2, LogLike.Mem2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem2$value + 6)
    }
    )
}
