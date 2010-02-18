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
    stop("Argument 1 of function SEM.MLE is not of class CLOngloneData")

  }else if( sum( subdata$N1+subdata$N2) == 0)
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
      sem2 = optim(theta.sem2,  LogLike.SEM2, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      if(verbose) print(sem2)
      hat.sem2 = exp(sem2$par)
      return(c(hat.sem2) )
    }
  )
}


##############################################
##    log likelihood for SEM1
##############################################
LogLike.SEM1 <- function(theta.sem1,  subdata) 
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
LogLike.SEM2 <- function(theta.sem2,  subdata) 
{
  ##theta.sem1=c(log(lambda1), log(lambda2), log(lambdae), log(sigma))
  lambda1 = exp(theta.sem2[1])
  lambda2 = exp(theta.sem2[2])
  lambdae = exp(theta.sem2[3])
  sigma   = exp(theta.sem2[4])
  lambda  = lambda1 + lambda2 + lambdae
  with(subdata,
     {
       return(
       sum( dnorm(N1, mean=lambda1/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) +
       sum( dnorm(N2, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) +
       sum( dnorm(Ne, mean=lambdae/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) )
              )
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
      sem1 = optim(theta.sem1,  LogLike.SEM1, method="BFGS", control=list(fnscale=-1),subdata=subdata)
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
      sem2 = optim(theta.sem2,  LogLike.SEM2, method="BFGS", control=list(fnscale=-1),subdata=subdata)
      return(-2*sem2$value + 8)
    }
  )
}
  

######################################################
##   Profile likelihood for lambda1 in MEM1
######################################################
Profile.SEM <- function(lambda1, subdata, SEM=1)
{
  ## lambda 1 is a vector 
  with(subdata,
       {
         return( sapply(lambda1, PartialLogLike.SEM, subdata=subdata, SEM=SEM))
      }
  )
}
######################################################
##   Partial likelihood for lambda1 in MEM1
##  optimized according other parameters
######################################################
PartialLogLike.SEM<- function(lambda1, subdata, SEM=1)
{
  with(subdata,
       {
         ##################
         ##CAS SEM 1
         if(SEM==1)
           {
             
             lambda.init  = 1/ mean(P)  * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             if(lambda1<lambda.init)
               {
                 lambda2.init = lambda.init - lambda1
               } else
             {
               lambda2.init=lambda1
               lambda.init=lambda2.init + lambda1  
             }
             sigma.init   = sqrt ( sum( (N1 - lambda1/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
               (N2+Ne - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 2 * length(N) ) )
             theta.sem1 = c( log(lambda2.init), log(sigma.init))
             ## definition of the function to be optimized only on lambda2 and p2
             partlog <- function(theta.sem1)
               {
                 return (LogLike.SEM1(c(log(lambda1), theta.sem1), subdata=subdata) )
               }
             optim.sem1 = optim(theta.sem1 , partlog, method="BFGS", control=list(fnscale=-1))
             return( optim.sem1$value )
           }else if(SEM==2)
             ####################################
             ## CAS SEM2
           {

             ## finding initial values lambda2, lambdae and sigma
             ## for optim
             ## if lambda1 
             lambda.init  = 1/ mean(P)  * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             if(lambda1<lambda.init)
               {
                 lambda2.init = sum(N2) /sum(N1+N2) *  ( lambda.init - lambda1)
                 lambdae.init = max(lambda.init - lambda1 -lambda2.init,0)
               } else
             {
               lambda2.init=lambda1/2
               lambdae.init=lambda1/2
               lambda.init=lambdae.init+lambda2.init+ lambda1
             }
             sigma.init   = sqrt (
               sum( (N1 - lambda1/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
               (N2 - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 +
               (Ne - lambdae.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 3 * length(N) ) )
             theta.sem2 = c( log(lambda2.init),log(lambdae.init), log(sigma.init))
             ## definition of the function to be optimized only on lambda2 and p2
             partlog <- function(theta.sem2)
               {
                 return (LogLike.SEM2(c(log(lambda1), theta.sem2), subdata=subdata) )
               }
             optim.sem2 = optim(theta.sem2 , partlog, method="BFGS", control=list(fnscale=-1))
             return( list(value=optim.sem2$value, lambda2=optim.sem2$par[2],lambdae=optim.sem2$par[2], sig=optim.sem2$par[3]))
           }
       }
  )
}


Prediction.SEM <-  function(theta, subdata)
{
  lambda1=theta[1]
  lambda2=theta[2]
  lambdae=theta[3]
  sigma=theta[4]
  lambda=lambda1+lambda2+lambdae
  with(subdata,
       {
         N1.pred= N*(1-exp(-lambda * P)) * lambda1 /lambda 
         N2.pred= N*(1-exp(-lambda * P)) * lambda2 /lambda
         Ne.pred= N*(1-exp(-lambda * P)) * lambdae /lambda
         return(list(N1.pred=N1.pred, N2.pred=N2.pred, Ne.pred=Ne.pred))
       }
       )
  }

