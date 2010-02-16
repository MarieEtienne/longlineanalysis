###########################################################################
## Method 1  : Non linear regression as proposed by Somerton and Kikkawa 
##
## BE CAUTIOUS : There is no meaning to compare AIC od SEM1 and SEM2 because SEM1 has less data(2 K) against 3k for SEM2 MEM2 MEM1
###########################################################################
## script to perform the Non linear regression of the exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## Y the dataframe which stored the data, longline in this code
##############################################################################

SEM.MLE <- function(subdata, sameN=T, sameP=TRUE, SEM=1, verbose=F) 
{
  ## initial values for minimization algorithm
  if(!is.CLonglineData(subdata))
  {
    stop("Argument 1 of function AnalysisMLE is not of class CLOngloneData")
  } else if( sum( subdata$N1+subdata$N2) == 0)
  {
    return(list( lambda = c( NA,NA,NA,NA ) ) ) ## no possible estimation
  } else if(SEM==1) 
  {
    return(Sem1MleGeneral(subdata, verbose=verbose))
  } else
  {
    return(Sem2MleGeneral(subdata, verbose=verbose))
  }
}

################################################
## Estimators for SEM1 general case
###############################################
Sem1MleGeneral <- function(subdata,verbose=F)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda2.init = sum( N2 + Ne) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda.init = lambda1.init + lambda2.init
      sigma.init   = sqrt ( sum( (N1 - lambda1.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
                    (N2+Ne - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 2 * length(N) ) )
      theta.sem1 = c( log(lambda1.init),log(lambda2.init),log(sigma.init))
      sem1 = optim(theta.sem1,  SEM1Loglike, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      if(verbose) print(sem1)
      hat.sem1 = exp(sem1$par)
      return(c(hat.sem1[1:2], NA, hat.sem1[3]) )
    }
  )
}
################################################
## Estimators for SEM2 general case
###############################################
Sem2MleGeneral <- function(subdata, verbose=F)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda2.init = sum( N2 ) / ( mean(P) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambdae.init = sum( Ne ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda.init  = lambda1.init + lambda2.init + lambdae.init
      sigma.init   = sqrt ( sum( (N1 - lambda1.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
                    (N2 - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 +
                    (Ne - lambdae.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 3 * length(N) ) )
      theta.sem2 = c( log(lambda1.init),log(lambda2.init),log(lambdae.init),log(sigma.init))
      sem2 = optim(theta.sem2,  SEM2Loglike, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      if(verbose) print(sem2)
      hat.sem2 = exp(sem2$par)
      return(c(hat.sem2) )
    }
  )
}


##############################################
##    log likelihood for SEM1
##############################################
SEM1Loglike <- function(theta.sem1,  subdata) 
{
  ##theta.sem1=c(log(lambda1), log(lambda2), log(sigma))
  lambda1 = exp(theta.sem1[1])
  lambda2 = exp(theta.sem1[2])
  sigma   = exp(theta.sem1[3])
  lambda  = lambda1 + lambda2
  with(subdata,
    {
      return( sum( dnorm(N1, mean=lambda1/lambda*N * (1-exp(-lambda * P) ), sd=sigma, log=T) )+
      sum( dnorm(N2+Ne, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) )
    }
  )
}
 
##############################################
##    log likelihood for SEM2
##############################################
SEM2Loglike <- function(theta.sem2,  subdata) 
{
  ##theta.sem1=c(log(lambda1), log(lambda2), log(lambdae), log(sigma))
  lambda1 = exp(theta.sem2[1])
  lambda2 = exp(theta.sem2[2])
  lambdae = exp(theta.sem2[3])
  sigma   = exp(theta.sem2[4])
  lambda  = lambda1 + lambda2 + lambdae
  with(subdata,
    {
      return( sum( dnorm(N1, mean=lambda1/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) )+
      sum( dnorm(N2, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) )   +
      sum( dnorm(Ne, mean=lambdae/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ))
    }
  )
}
 

######################################################
##   SEM1 AIC
######################################################
AIC.Sem1 <- function(subdata)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda2.init = sum( N2 + Ne) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda.init  = lambda1.init + lambda2.init
      sigma.init   = sqrt ( sum( (N1 - lambda1.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
                    (N2+Ne - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 2 * length(N) ) )
      theta.sem1 = c( log(lambda1.init),log(lambda2.init),log(sigma.init))
      sem1 = optim(theta.sem1,  SEM1Loglike, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      return(-2*sem1$value+6)
    }
  )
}
  

######################################################
##   SEM2 AIC
######################################################
AIC.Sem2 <- function(subdata)
{
  with(subdata,
    {
      lambda1.init = sum( N1 ) / ( mean(P) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda2.init = sum( N2 ) / ( mean(P) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambdae.init = sum( Ne ) / ( mean(P) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      lambda.init  = lambda1.init + lambda2.init + lambdae.init
      sigma.init   = sqrt ( sum( (N1 - lambda1.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
                    (N2 - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 +
                    (Ne - lambdae.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 3 * length(N) ) )
      theta.sem2 = c( log(lambda1.init),log(lambda2.init),log(lambdae.init),log(sigma.init))
      sem2 = optim(theta.sem2,  SEM2Loglike, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      return(-2*sem2$value + 8)
    }
  )
}
  

