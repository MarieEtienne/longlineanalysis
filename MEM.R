##############################################################
## Method 2  : ML estimation 
##############################################################
## script to perform the ML estimation of the complete exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## subdata the dataframe which stored the data to be analyzed
##############################################################################

MEM.MLE <- function(subdata, tune1=3, sameP=TRUE, MEM=1) 
{
  ## initial values for minimization algorithm
  if(!is.CLonglineData(subdata))
  {
    stop("Argument 1 of function MEM.MLE is not of class CLOngloneData")
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
      optim.mem1 = optim(theta.mem1, LogLike.MEM1, method="BFGS", control=list(fnscale=-1),  subdata=subdata)
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
      
      optim.mem2 = optim(theta.mem2, LogLike.MEM2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
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
##   MEM LogLikelihood
######################################################
LogLike.MEM1 <- function(theta.mem1,  subdata)
{
  ## theta.mem1 = c(log(lambda1), log(lambda2), logit(p2))
  return( LogLike( c(theta.mem1[1:2], -Inf, theta.mem1[3]), subdata=subdata ) )
}

LogLike.MEM2 <- function(theta.mem2,  subdata)
{
  ## theta.mem2 = c(log(lambda1), log(lambda2), logit(p))
  return( LogLike( c(theta.mem2[1:2], theta.mem2[3], theta.mem2[3] ), subdata=subdata ) )
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
      optim.mem1 = optim(theta.mem1, LogLike.MEM1, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem1$value + 6)
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
      
      optim.mem2 = optim(theta.mem2, LogLike.MEM2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem2$value + 6)
    }
    )
}


######################################################
##   Profile likelihood for lambda1 in MEM1
######################################################
Profile.MEM <- function(lambda1, subdata, MEM=1)
{
  ## lambda 1 is a vector 
  with(subdata,
       {
         return( sapply(lambda1, PartialLogLike, subdata=subdata, MEM=MEM))
      }
  )
}
######################################################
##   Partial likelihood for lambda1 in MEM1
##  optimized according other parameters
######################################################
PartialLogLike.MEM<- function(lambda1, subdata, MEM=1)
{
  with(subdata,
       {
         if(MEM==1)
           {
             lambda2.init = sum( N2 + Ne) / ( mean(P) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             p2           = sum( Ne )/ sum( Ne + N2 )
             theta.mem1 = c( log(lambda2.init), log(p2/(1-p2)))
             ## definition of the function to be optimized only on lambda2 and p2
             partlog <- function(theta.mem1)
               {
                 return (LogLike.MEM1(c(log(lambda1), theta.mem1), subdata=subdata) )
               }
             optim.mem1 = optim(theta.mem1 , partlog, method="SANN", control=list(fnscale=-1))
             return( optim.mem1$value )
           }else if(MEM==2)
           {
             lambda2.init = sum( N2 ) / ( mean(P) * sum( (N1 + N2) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             p            = sum( Ne )/ sum( Ne + N2  + N1)
             theta.mem2   = c( log(lambda2.init), log(p/(1-p)))
             ## definition of the function to be optimized only on lambda2 and p2
             partlog <- function(theta.mem2)
               {
                 return (LogLike.MEM2(c(log(lambda1), theta.mem2), subdata=subdata) )
               }
             optim.mem2 = optim(theta.mem2 , partlog, method="BFGS", control=list(fnscale=-1))
             return( optim.mem2$value )
           }
       }
  )
}


#################################################
##Prediction of expected values
##
#################################################
Prediction.MEM <-  function(theta, subdata)
{
  lambda1=theta[1]
  lambda2=theta[2]
  lambda=lambda1+lambda2
  p1=theta[3]
  p2=theta[4]
  with(subdata,
       {
         N1.pred= N*(1-exp(-lambda * P)) * lambda1 /lambda *(1-p1)
         N2.pred= N*(1-exp(-lambda * P)) * lambda2 /lambda *(1-p2)
         Ne.pred= N*(1-exp(-lambda * P)) * (lambda1*p1 +lambda2*p2) /lambda
         return(list(N1.pred=N1.pred, N2.pred=N2.pred, Ne.pred=Ne.pred))
       }
       )
  }


################################################
## Covariance matrix asymptotic
## only valid if all P are similar
###############################################

MEM.Cov <- function(subdata, MEM=1)
  {
    P=mean(subdata$P)
    theta=MEM.MLE(subdata, MEM=MEM)
    lambda1   = theta[1]
    lambda2   = theta[2]
    p1        = theta[3]
    p2        = theta[4]
    lambda    = lambda1+lambda2
    moinsexp = 1-exp(-lambda*P)

    if(MEM==1)
      {
        mat= matrix(c(lambda1*lambda2 / moinsexp + moinsexp *lambda1^2 / ( P^2 * exp(-lambda *P) * lambda^2),
          -lambda1*lambda2 / moinsexp + moinsexp *lambda1*lambda2 / ( P^2 * exp(-lambda *P) * lambda^2),
          0,
          -lambda1*lambda2 / moinsexp + moinsexp *lambda1*lambda2 / ( P^2 * exp(-lambda *P) * lambda^2),
          lambda1*lambda2 / moinsexp + moinsexp *lambda2^2 / ( P^2 * exp(-lambda *P) * lambda^2),
          0,
          0,
          0,
          lambda*p2*(1-p2)/ (lambda2 * moinsexp)), byrow=T, ncol=3)
      }else if( MEM==2)
      {
      mat= matrix(c(lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1^2 / ( P^2 * exp(-lambda *P) * lambda^2),
          -lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1*lambda2 / ( P^2 * exp(-lambda *P) * lambda^2),
          0,
          -lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1*lambda2 / ( P^2 * exp(-lambda *P) * lambda^2),
          lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda2^2 / ( P^2 * exp(-lambda *P) * lambda^2),
          0,
          0,
          0,
          p2*(1-p2)/ (moinsexp)), byrow=T, ncol=3)
        
      }
    
    
    return(1/sum(subdata$N) * mat)
  }
