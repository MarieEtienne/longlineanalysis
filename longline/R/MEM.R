##############################################################
## Method 2  : ML estimation 
##############################################################

## script to perform the ML estimation of the complete exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## subdata the dataframe which stored the data to be analyzed
##############################################################################


# MEM.MLE
#' MEM.MLE computes MLE for data with several years and area, one lambda1, lambda2, pe computes for each(year,area)
#' 
#' @param llData a object of class longline
#' @param MEM the type of model (possible values 1 or 2)
#' @param sameS, should all soaktime be considered as equal ?
#' 
#' @export
#' @return a vetor estimates of hat.lmabdaT, hat.lambda.NT, hat.pT, hat.pNT
#
MEM.MLE <- function(llData, MEM=1, sameS=T) 
{
  adjustFact <- with(llData, exists("Fact2"))
  if(sameS){
    llData$sameS=T
   }
  if(adjustFact) {
    llData$Fact1 <- as.factor(with(llData, paste0(Fact1, "-", Fact2 )))
    llData$Fact2 <- NULL
    llData$NFact2 <- NULL
  }      

  res<- lapply(levels(llData$Fact1), function(d){
    subdata <- ExtractLongline(llData, fact1 = d)
    tt <-MEM.MLE.YA(subdata=subdata, MEM=MEM)
    return(tt)
  })
  
  estimate <- do.call("rbind", res)
  colnames(estimate) <- c("hat.lambda1", "hat.lambda2", "hat.pT", "hat.pNT")
  estimate <- as.data.frame(estimate)
  estimate$Group <- levels(llData$Fact1)
  return(estimate)
}




#### mem.mle.YA computes MEM MLE for one year in one area
MEM.MLE.YA <- function(subdata, MEM=1) 
{
  ## initial values for minimization algorithm
  if(!is.longline(subdata))
  {
    stop("Argument 1 of function MEM.MLE is not of class longline")
  } else if( sum( subdata$NT+subdata$NNT) == 0)
  {
    return(list( lambda = c( , NA,NA, NA,NA ) ) ) ## no possible estimation
  } else if(subdata$sameS)     
  {
    ne.missing <- (sum(is.na(subdata$Ne))==subdata$NData)
    if(MEM==1) 
    {
      return(Mem1MleSameP(subdata, ne.missing=ne.missing))
    }else if(MEM==2)
    {
      return(Mem2MleSameP(subdata, ne.missing=ne.missing))
    }             
  } else
  {
    ne.missing <- (sum(is.na(subdata$Ne))==subdata$NData)
    if(MEM==1) 
    {
      return(Mem1MleGeneral(subdata, ne.missing=ne.missing))
    } else
    {
      return(Mem1MleGeneral(subdata, ne.missing=ne.missing))
    }
  }
}


################################################
## Estimators for MEM1 with shared soaktime P
###############################################
Mem1MleSameP <- function(subdata, ne.missing=F)
{
  with(subdata,
    {
      hat.lambda1=sum( NT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.lambda1 <- max(1e-7, hat.lambda1)
      if(ne.missing)
      {
        hat.lambda2=sum( NNT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
        hat.p2 = 0.5
      }else
       {
         hat.lambda2=sum( NNT + Ne) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
         hat.p2 = sum(Ne)/ sum(Ne+NNT)
         
       } 
      
      hat.p1 =0
      return(c(hat.lambda1, hat.lambda2, hat.p1, hat.p2))
    }
  )
}

################################################
## Estimators for MEM2 with shared soaktime P
###############################################
Mem2MleSameP <- function(subdata, ne.missing=F)
{
  with(subdata,
    {
      hat.lambda1=sum( NT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      hat.lambda2=sum( NNT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
     if(ne.missing)
     {
       hat.p2=0.5
       hat.p1=0.5
     }else
     {
        hat.p2 = sum( Ne )/ sum( N - Nb )
        hat.p1 =  hat.p2 
     }
      return(c(hat.lambda1, hat.lambda2, hat.p1, hat.p2))
    }
  )
}

################################################
## Estimators for MEM1 general case
###############################################
Mem1MleGeneral <- function(subdata, ne.missing=F)
{
  with(subdata,
    {                                    
      lambda1.init = sum( NT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb)) 
      if(ne.missing)
      {
      lambda2.init = sum( NNT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = 0.5
      }else
      {
        lambda2.init = sum( NNT + Ne ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
        p            = sum( Ne )/ sum( Ne + NNT)
      }
      theta.mem1 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p) ) )
      if(ne.missing)
        par.mem1   = theta.mem1
      else
        {
        optim.mem1 = optim(theta.mem1, LogLike.MEM1, method="BFGS", control=list(fnscale=-1),  subdata=subdata)
        par.mem1   = optim.mem1$par 
        }
      return( c(exp(par.mem1[1:2]), 0, exp(par.mem1[3]) / (1 + exp(par.mem1[3]) )) )
    }
  )
}


################################################
## Estimators for MEM2 general case
###############################################
Mem2MleGeneral <- function(subdata, ne.missing=F)
{
  with(subdata,
    {
      lambda1.init = sum( NT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( NNT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      if(ne.missing)
      {
        p            = 0.5
      }else
      {
        p            = sum( Ne )/ sum( Ne + NNT + NT )
      }
      theta.mem2 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      if(ne.missing)
        par.mem2   = theta.mem2
      else
        {
          optim.mem2 = optim(theta.mem2, LogLike.MEM2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
        par.mem2   = optim.mem2$par 
        }
      return( c(par.mem2, par.mem2[3]) )
    }
  )
}



######################################################
##   General LogLikelihood 
######################################################
## only when Ne is not missing
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
        sum( lgamma(N+1)-lgamma(Nb+1)-lgamma(NT+1)-lgamma(NNT+1)-lgamma(Ne+1) ) - sum(lambda * S * Nb) +
        sum( (N-Nb) * log (1 - exp( -lambda * S) ) ) + sum(NT * log(  lambda1 / lambda * ( 1-p1) ) )+ 
        sum( NNT * log(  lambda2 / lambda * ( 1-p2) ) ) + sum ( Ne * log(  (lambda1 *p1 + lambda2 * p2) / lambda ) )
        )
    }
  )
}

######################################################
##   MEM LogLikelihood
######################################################
## only when Ne is not missing
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
## only when Ne is not missing
AIC.Mem1 <- function(subdata)
{
   with(subdata,
    {
      lambda1.init = sum( NT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( NNT ) / ( mean(S) * sum( (N - Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + NNT )
      theta.mem1 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      optim.mem1 = optim(theta.mem1, LogLike.MEM1, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem1$value + 6)
    }
    )
}

######################################################
##   MEM2 AIC
######################################################
## only when Ne is not missing
AIC.Mem2 <- function(subdata)
{
   with(subdata,
    {
      lambda1.init = sum( NT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb)) 
      lambda2.init = sum( NNT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) since all p are supposed to be equal or can be considered as equal
      p            = sum( Ne )/ sum( Ne + NNT + NT )
      theta.mem2 = c(log(lambda1.init), log(lambda2.init), log(p/(1-p)))
      
      optim.mem2 = optim(theta.mem2, LogLike.MEM2, method="BFGS", control=list(fnscale=-1), subdata=subdata)
      return( -2 * optim.mem2$value + 6)
    }
    )
}


######################################################
##   Profile likelihood for lambda1 in MEM1
######################################################
## only when Ne is not missing
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
## only when Ne is not missing
PartialLogLike.MEM<- function(lambda1, subdata, MEM=1)
{
  with(subdata,
       {
         if(MEM==1)
           {
             lambda2.init = sum( NNT + Ne) / ( mean(S) * sum( (N- Nb) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             p2           = sum( Ne )/ sum( Ne + NNT )
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
             lambda2.init = sum( NNT ) / ( mean(S) * sum( (NT + NNT) ) ) * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
             p            = sum( Ne )/ sum( Ne + NNT  + NT)
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
## only when Ne is not missing
Prediction.MEM <-  function(theta, subdata)
{
  lambda1=theta[1]
  lambda2=theta[2]
  lambda=lambda1+lambda2
  p1=theta[3]
  p2=theta[4]
  with(subdata,
       {
         NT.pred= N*(1-exp(-lambda * P)) * lambda1 /lambda *(1-p1)
         NNT.pred= N*(1-exp(-lambda * P)) * lambda2 /lambda *(1-p2)
         Ne.pred= N*(1-exp(-lambda * P)) * (lambda1*p1 +lambda2*p2) /lambda
         return(list(NT.pred=NT.pred, NNT.pred=NNT.pred, Ne.pred=Ne.pred))
       }
       )
  }


################################################
## Covariance matrix asymptotic
## only valid if all P are similar
###############################################
## only when Ne is not missing
MEM.Cov <- function(subdata, MEM=1)
  {
    S=mean(subdata$S)
    theta=MEM.MLE(subdata, MEM=MEM)
    lambda1   = theta[1]
    lambda2   = theta[2]
    p1        = theta[3]
    p2        = theta[4]
    lambda    = lambda1+lambda2
    moinsexp = 1-exp(-lambda*S)

    if(MEM==1)
      {
        mat= matrix(c(lambda1*lambda2 / moinsexp + moinsexp *lambda1^2 / ( S^2 * exp(-lambda *S) * lambda^2),
          -lambda1*lambda2 / moinsexp + moinsexp *lambda1*lambda2 / ( S^2 * exp(-lambda *S) * lambda^2),
          0,
          -lambda1*lambda2 / moinsexp + moinsexp *lambda1*lambda2 / ( S^2 * exp(-lambda *S) * lambda^2),
          lambda1*lambda2 / moinsexp + moinsexp *lambda2^2 / ( S^2 * exp(-lambda *S) * lambda^2),
          0,
          0,
          0,
          lambda*p2*(1-p2)/ (lambda2 * moinsexp)), byrow=T, ncol=3)
      }else if( MEM==2)
      {
      mat= matrix(c(lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1^2 / ( S^2 * exp(-lambda *S) * lambda^2),
          -lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1*lambda2 / ( S^2 * exp(-lambda *S) * lambda^2),
          0,
          -lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda1*lambda2 / ( S^2 * exp(-lambda *S) * lambda^2),
          lambda1*lambda2 / ((1-p2)* moinsexp) + moinsexp *lambda2^2 / ( S^2 * exp(-lambda *S) * lambda^2),
          0,
          0,
          0,
          p2*(1-p2)/ (moinsexp)), byrow=T, ncol=3)
        
      }
      
    
    return(1/sum(subdata$N) * mat)
  }


  MEM.Fisher <- function(theta, S, MEM=1)
  {
    lambda1   = theta[1]
    lambda2   = theta[2]
    p1        = theta[3]
    p2        = theta[4] 
    p         = p1               
    lambda    = lambda1+lambda2
    moinsexp = 1-exp(-lambda*S)

    if(MEM==1)
      {
        mat= matrix(c(S*S * exp(-lambda * S) / moinsexp + moinsexp /( lambda) *(1/lambda1 - 1/lambda),
                      S*S * exp(-lambda * S) / moinsexp + moinsexp /( lambda*lambda),
                      0,
                      S*S * exp(-lambda * S) / moinsexp + moinsexp /( lambda*lambda),
                      S*S * exp(-lambda * S) / moinsexp + moinsexp /( lambda) *(1/lambda2 - 1/lambda),
                      0,
                      0,
                      0,
                      (lambda2 * moinsexp)/(lambda  *p2*(1-p2) ))
        , byrow=T, ncol=3)
      }else if( MEM==2)
      {
        mat= matrix(c(S*S * exp(-lambda * S) / moinsexp + (1-p) * moinsexp /( lambda * lambda ) *(lambda2 /lambda1),
                      S*S * exp(-lambda * S) / moinsexp - (1-p)  * moinsexp /( lambda*lambda),
                      0,
                      S*S * exp(-lambda * S) / moinsexp - (1-p)  * moinsexp /( lambda*lambda),
                      S*S * exp(-lambda * S) / moinsexp + (1-p) * moinsexp /( lambda * lambda ) *(lambda1 /lambda2),
                      0,
                      0,
                      0,
                      (moinsexp)/(p*(1-p) ))
        , byrow=T, ncol=3) 
      }
    return(mat)
  }
