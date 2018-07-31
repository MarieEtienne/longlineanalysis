## 7/09/2016

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

# SEM.MLE
#' SEM.MLE computes MLE for data with several years and area, one lambda1, lambda2, computed for each(year,area)
#' 
#' @param llData a object of class longline
#' @param sameS, should all soaktime be considered as equal ?
#' @param sameN, should all initial number of baits be considered as equal ?
#' @param sameSigma, should the variance be chosen equal for target and Non Target
#' @param verbose, for debugging purpose
#' @param SEM, the version of the SEM model to be calculated : possible value 1 (if NE and NNT are grouped together, same variance), 2 (NE and NNT split into 2 categories, Same variance), 3(NE and NNT together, two different variance sigmaT and SigmaNT)
#' 
#' @export
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT

SEM.MLE <- function(llData, sameN=TRUE, sameS=TRUE, SameSigma=T, verbose=F, SEM=1) 
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
    tt <-SEM.MLE.YA(subdata=subdata, sameN=TRUE, sameP=TRUE, SameSigma=T, verbose=F, SEM=SEM)
    return(tt)
  })
  estimate <- do.call("rbind", res)
  colnames(estimate) <- c("hat.lambda1", "hat.lambda2", "hat.sigma1", "hat.sigma2")
  estimate <- as.data.frame(estimate)
  estimate$Group <- levels(llData$Fact1)
  return(estimate)
  
}


# SEM.MLE.YA
#' SEM.MLE.YA computes MLE for data for one specific year, one lambda1, lambda2, computed for each(area)
#' 
#' @param llData a object of class longline
#' @param sameS, should all soaktime be considered as equal ?
#' @param sameN, should all initial number of baits be considered as equal ?
#' @param sameSigma, should the variance be chosen equal for target and Non Target
#' @param verbose, for debugging purpose
#' @param SEM, the version of the SEM model to be calculated : possible value 1 (if NE and NNT are grouped together, same variance), 2 (NE and NNT split into 2 categories, Same variance), 3(NE and NNT together, two different variance sigmaT and SigmaNT)
#' 
#' @export
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT

SEM.MLE.YA <- function(subdata, sameN=TRUE, sameS=TRUE, SameSigma=1, verbose=F, SEM=SEM) {
  ## initial values for minimization algorithm
  if(!is.longline(subdata))
  {
    stop("Argument 1 of function SEM.MLE is not of class longline")
  } else if( sum( subdata$NT+subdata$NNT) == 0)
  {
    return(list( theta = c( NA,NA, NA,NA,NA ) ) ) ## no possible estimation
  } else if(subdata$sameS & (SEM==1|SEM==2) )     
  {
    if(sameSigma & SEM==1) 
    {
      return(Sem1MleSameS(subdata))
    }else if(sameSigma & SEM==2) 
    {
      return(Sem2MleSameS(subdata))
    }
  } else
  {
    if(sameSigma & SEM==1) 
    {
      return(Sem1MleGeneral(subdata))
    } else if (sameSigma & SEM==2)
    {
      return(Sem2MleGeneral(subdata))
    } else
    {
      return(Sem3MleGeneral(subdata))
    }
  }
}


# Sem1MleSameS
#' Sem1MleSameS computes MLE for data for one specific year, and one specific area when all longlinesets share the same soaktime
#' if the soaktime are different a mean soaktime is computed
#' 
#' @param data a object of class longline
#' 
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT

Sem1MleSameS <- function(data){
  with(subdata, 
       {
         P.m=mean(P)
         N.m =mean(N)
         lambda1 = sum(NT) / sum (N-Nb) * 1/P.m *log(sum(N) /sum(Nb) )
         lambda2 = sum(NNT+Ne) / sum (N-Nb) * 1/P.m *log(sum(N) /sum(Nb) )
         lambda=lambda1+lambda2
         sigma = sqrt(1/(2*length(N))*sum ( (NT-N*lambda1 /lambda *(1-exp(-lambda *P)))^2 + (NNT+Ne-N*lambda2 /lambda *(1-exp(-lambda *P)))^2)  )
         return(c( lambda1,lambda2,NA,sigma,NA )  ) ## no possible estimation 
       })
}


# Sem2MleSameS
#' Sem2MleSameS computes MLE for data for one specific year, and one specific area when all longlinesets share the same soaktime
#' if the soaktime are different a mean soaktime is computed
#' 
#' @param data a object of class longline
#' 
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT

Sem2MleSameS <- function(data){
  with(subdata, 
       {
         P.m=mean(P)
         N.m =mean(N)
         lambda1 = sum(NT) / sum (N-Nb) * 1/P.m *log(sum(N) /sum(Nb) )
         lambda2 = sum(NNT) / sum (N-Nb) * 1/P.m *log(sum(N) /sum(Nb) )
         lambdae = sum(Ne) / sum (N-Nb) * 1/P.m *log(sum(N) /sum(Nb) )
         lambda=lambda1+lambda2+lambdae
         sigma = sqrt(1/(2*length(N))*sum ( (NT-N*lambda1 /lambda *(1-exp(-lambda *P)))^2 + 
                                              (NNT-N*lambda2 /lambda *(1-exp(-lambda *P)))^2 +
                        (Ne-N*lambdae /lambda *(1-exp(-lambda *P)))^2)  )
         return(c( lambda1,lambda2,NA,sigma,NA )  ) ## no possible estimation 
       })
}





# Sem1MleGeneral
#' Sem1MleSameS computes MLE for data for one specific year, and one specific area 
#' 
#' @param data a object of class longline
#' @param verbose for debugging purpose
#' 
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT
Sem1MleGeneral <- function(subdata,verbose=F)
{
  with(subdata,
       {      
         theta.init <- Sem1MleSameS(data)
         theta.sem1 = c( log(theta.init[1]),log(theta.init[2]), log(theta.init[4]))
         sem1 = optim(theta.sem1,  LogLike.SEM1, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         if(verbose) print(sem1)
         hat.sem1 = exp(sem1$par)
         return(c(hat.sem1[1:2], NA, hat.sem1[3],NA) )
       }
  )
}


## Sem2MleGeneral
#' Sem2MleSameS computes MLE for data for one specific year, and one specific area 
#' 
#' @param data a object of class longline
#' @param verbose for debugging purpose
#' 
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT
Sem2MleGeneral <- function(subdata, verbose=F)
{
  with(subdata,
       {
         theta.init <- Sem2MleSameS(data)
         theta.sem2 = c( log(theta.init[1]),log(theta.init[2]), log(theta.init[3]),log(theta.init[4]))
         sem2 = optim(theta.sem2,  LogLike.SEM2, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         if(verbose) print(sem2)
         hat.sem2 = exp(sem2$par)
         return(c(hat.sem2[1:4], NA) )
       }
  )
}


# Sem3MleGeneral
#' Sem3MleSameS computes MLE for data for one specific year, and one specific area 
#' 
#' @param data a object of class longline
#' @param verbose for debugging purpose
#' 
#' @return a vetor estimates of hat.lambdaT, hat.lambda.NT,  hat.lambdaNe, hat.sigmaT, hat.sigmaNT
Sem3MleGeneral <- function(subdata, verbose=F)
{
  with(subdata,
       {
         theta.init <- Sem1MleSameS(data)
         theta.sem3 = c( log(theta.init[1]),log(theta.init[2]),log(theta.init[4]), log(theta.init[5]))
         sem3 = optim(theta.sem3,  LogLike.SEM3, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         if(verbose) print(sem3)
         hat.sem3 = exp(sem3$par)
         return(c(hat.sem2[1:2], NA, hat.sem2[3:4]) )
       }
  )
}

# LogLike.SEM1 
#' LogLike.SEM1 the log likelihood given data for the parameter theta
#' 
#' @param theta.sem1 parameter relevant for SEM 1 model in log scale
#' @param data a object of class longline
#' 
#' @return the log likelihood for SEM1
# LogLike.SEM1 <- function(theta.sem1,  subdata) 
# {
#   ##theta.sem1=c(log(lambda1), log(lambda2), log(sigma))
#   lambda1 = exp(theta.sem1[1])
#   lambda2 = exp(theta.sem1[2])
#   sigma   = exp(theta.sem1[3])
#   lambda  = lambda1 + lambda2
#   with(subdata,
#        {
#          return( sum( dnorm(NT, mean=lambda1/lambda*N * (1-exp(-lambda * P) ), sd=sigma, log=T) )+
#                    sum( dnorm(NNT+Ne, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) )
#        }
#   )
# }
# 
# LogLike.SEM2 
#' LogLike.SEM2 the log likelihood given data for the parameter theta
#' 
#' @param theta.sem2 parameter relevant for SEM 2 model in log scale
#' @param data a object of class longline
#' 
#' @return the log likelihood for SEM2
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
           sum( dnorm(NT, mean=lambda1/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) +
             sum( dnorm(NNT, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) ) +
             sum( dnorm(Ne, mean=lambdae/lambda*N * (1-exp(-lambda * P)), sd=sigma, log=T) )
         )
       }
  )
 }
 
 
 
 # LogLike.SEM3 
 #' LogLike.SEM3 the log likelihood given data for the parameter theta
 #' 
 #' @param theta.sem3 parameter relevant for SEM 3 model in log scale
 #' @param data a object of class longline
 #' 
 #' @return the log likelihood for SEM3
 LogLike.SEM3 <- function(theta.sem3,  subdata) 
{
  ##theta.sem1=c(log(lambda1), log(lambda2), log(lambdae), log(sigma))
  lambda1 <-  exp(theta.sem3[1])
  lambda2 <-  exp(theta.sem3[2]) ## theta.Sem3 contains lambdaempty whic doesn't make sense here
  sigma1  <-  exp(theta.sem3[4])
  sigma2  <-  exp(theta.sem3[5])
  
  lambda  = lambda1 + lambda2 
  with(subdata,
       {
         return(
           sum( dnorm(NT, mean=lambda1/lambda*N * (1-exp(-lambda * P)), sd=sigma1, log=T) ) +
             sum( dnorm(NNT, mean=lambda2/lambda*N * (1-exp(-lambda * P)), sd=sigma2, log=T) )
         )
       }
  )
 }
 
 
 # AIC.SEM
 #' AIC.SEM computes AIC for SEM Model
 #' 
 #' @param data the longline data
 #' @param SEM the SEM version
 #' 
 #' @return the AIC criterion
 
 AIC.SEM <- function(subdata, SEM=SEM)
{
  with(subdata,
       {if(SEM==1){
         theta.init <- Sem1MleSameS(data)
         theta.sem1 = c( log(theta.init[1]),log(theta.init[2]), log(theta.init[4]))
         sem1 = optim(theta.sem1,  LogLike.SEM1, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         return(-2*sem1$value+6)
       } else if (SEM==2){
         theta.init <- Sem2MleSameS(data)
         theta.sem2 = c( log(theta.init[1]),log(theta.init[2]), log(theta.init[3]),log(theta.init[4]))
         sem2 = optim(theta.sem2,  LogLike.SEM2, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         return(-2*sem2$value+8)
       } else if(SEM==3){
         theta.init <- Sem1MleSameS(data)
         theta.sem3 = c( log(theta.init[1]),log(theta.init[2]),log(theta.init[4]), log(theta.init[5]))
         sem3 = optim(theta.sem3,  LogLike.SEM3, method="BFGS", control=list(fnscale=-1),subdata=subdata)
         return(-2*sem2$value+8)
       } else{
         cat('SEM Version should be 1, 2 or 3\n')
         return(NA)
       }
       }
  )
}
# # ######################################################
# # ##   Profile likelihood for lambda1 in MEM1
# # ######################################################
# # Profile.SEM <- function(lambda1, subdata, SEM=1)
# # {
# #   ## lambda 1 is a vector 
# #   with(subdata,
# #        {
# #          return( sapply(lambda1, PartialLogLike.SEM, subdata=subdata, SEM=SEM))
# #        }
# #   )
# # }
# # ######################################################
# # ##   Partial likelihood for lambda1 in MEM1
# # ##  optimized according other parameters
# # ######################################################
# # PartialLogLike.SEM<- function(lambda1, subdata, SEM=1)
# # {
# #   with(subdata,
# #        {
# #          ##################
# #          ##CAS SEM 1
# #          if(SEM==1)
# #          {
# #            
# #            lambda.init  = 1/ mean(P)  * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
# #            if(lambda1<lambda.init)
# #            {
# #              lambda2.init = lambda.init - lambda1
# #            } else
# #            {
# #              lambda2.init=lambda1
# #              lambda.init=lambda2.init + lambda1  
# #            }
# #            sigma.init   = sqrt ( sum( (NT - lambda1/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
# #                                         (NNT+Ne - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 2 * length(N) ) )
# #            theta.sem1 = c( log(lambda2.init), log(sigma.init))
# #            ## definition of the function to be optimized only on lambda2 and p2
# #            partlog <- function(theta.sem1)
# #            {
# #              return (LogLike.SEM1(c(log(lambda1), theta.sem1), subdata=subdata) )
# #            }
# #            optim.sem1 = optim(theta.sem1 , partlog, method="BFGS", control=list(fnscale=-1))
# #            return( optim.sem1$value )
# #          }else if(SEM==2)
# #            ####################################
# #          ## CAS SEM2
# #          {
# #            
# #            ## finding initial values lambda2, lambdae and sigma
# #            ## for optim
# #            ## if lambda1 
# #            lambda.init  = 1/ mean(P)  * log( sum(N) / sum(Nb))  # P=mean(P) : initial guess
# #            if(lambda1<lambda.init)
# #            {
# #              lambda2.init = sum(NNT) /sum(NT+NNT) *  ( lambda.init - lambda1)
# #              lambdae.init = max(lambda.init - lambda1 -lambda2.init,0)
# #            } else
# #            {
# #              lambda2.init=lambda1/2
# #              lambdae.init=lambda1/2
# #              lambda.init=lambdae.init+lambda2.init+ lambda1
# #            }
# #            sigma.init   = sqrt (
# #              sum( (NT - lambda1/lambda.init*N *(1-exp(-lambda.init*P)))^2 + 
# #                     (NNT - lambda2.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 +
# #                     (Ne - lambdae.init/lambda.init*N *(1-exp(-lambda.init*P)))^2 ) / ( 3 * length(N) ) )
# #            theta.sem2 = c( log(lambda2.init),log(lambdae.init), log(sigma.init))
# #            ## definition of the function to be optimized only on lambda2 and p2
# #            partlog <- function(theta.sem2)
# #            {
# #              return (LogLike.SEM2(c(log(lambda1), theta.sem2), subdata=subdata) )
# #            }
# #            optim.sem2 = optim(theta.sem2 , partlog, method="BFGS", control=list(fnscale=-1))
# #            return( list(value=optim.sem2$value, lambda2=optim.sem2$par[2],lambdae=optim.sem2$par[2], sig=optim.sem2$par[3]))
# #          }
# #        }
# #   )
# # }
# # 
# # 
# # Prediction.SEM <-  function(theta, subdata)
# # {
# #   lambda1=theta[1]
# #   lambda2=theta[2]
# #   lambdae=theta[3]
# #   sigma=theta[4]
# #   lambda=lambda1+lambda2+lambdae
# #   with(subdata,
# #        {
# #          NT.pred= N*(1-exp(-lambda * P)) * lambda1 /lambda 
# #          NNT.pred= N*(1-exp(-lambda * P)) * lambda2 /lambda
# #          Ne.pred= N*(1-exp(-lambda * P)) * lambdae /lambda
# #          return(list(NT.pred=NT.pred, NNT.pred=NNT.pred, Ne.pred=Ne.pred))
# #        }
# #   )
# # }
# # 
