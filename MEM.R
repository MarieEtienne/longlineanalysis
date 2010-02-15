##############################################################
## Method 2  : ML estimation 
##############################################################
## script to perform the ML estimation of the complete exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## data the dataframe which stored the data to be analyzed
## tune1 a parameter to tune due to the problem with numerical 
##       minization, one may try value as 3,4, 5, 2
##############################################################################




## function to minimize in order to obtain the MLestimator for the exponential model, 
## with several longlines, with different soaktimes
EquationNormaleLambda1 <- function(lambda2, Y)
{
  cn <- 1 + with(Y, sum(N1) / sum(N2) )
  res <- with ( Y, - sum (P * Nb)  + 
                           sum(  (N- Nb) * P * exp( - lambda2  * cn *P ) / 
                                 ( 1 - exp( - lambda2  * cn *P ) ) 
                               )
              )
    return( log(res^2+1) )
}


AnalysisMLE <- function(subdata, tune1=3) 
{
  ## initial values for minimization algorithm
  if(!is.CLonglineData(subdata))
    {
      stop("Argument 1 of function AnalysisMLE is not of class CLOngloneData")
    } else
    {
      if( sum( subdata$N1+subdata$N2 ) == 0)
        {
          return(list( lambda = c( NA,NA ) ) )
        } else
        {
          
          p <- mean(subdata$P)
          algo.min.init <- with(subdata,sum(N1)/sum(N1+N2) * log( sum(N) /sum(Nb) )/p)
          value.cur <- EquationNormaleLambda1(algo.min.init, subdata)
          algo.min <-  nlminb(algo.min.init, EquationNormaleLambda1, lower = 0.0000001, upper = 0.01, Y=subdata)
          iter=0
          
          ##print( paste(" Objective=", algo.min$objective, ", Par=", algo.min$par, sep="")) 
          
          while(iter < 10 & algo.min$objective< value.cur)
            {
              iter <- iter + 1
              algo.min.init<- algo.min$par
              value.cur <- algo.min$objective
              algo.min <-  nlminb(algo.min.init, EquationNormaleLambda1, lower = algo.min.init - tune1 * algo.min.init / 10 , upper =algo.min.init +  tune1 * algo.min.init/10, Y=subdata) 
            }
          lambda2.est <- algo.min$par
          lambda1.est <- with(subdata, sum(N1) / sum (N2) ) * lambda2.est
          ##print(paste("lambda1.est=",lambda1.est,sep=""))
          ##print(paste("lambda2.est=", lambda2.est,sep=""))
          return(list(lambda = c(lambda1.est, lambda2.est)) )
        }
    }
}


#################################################################
##
## Code for the log likelihood profile, for confidence interval
## construction
##################################################################

MEMLoglike <- function(lambda2,lambda1,Y)
  {
    lambda <- lambda1 + lambda2
    l <- with(Y,
              sum(- lambda * Nb * P  +
                  ( N1 + N2 + Ne) * log(1- exp(-lambda * P)) ) +
              sum(N1) * log(lambda1) + sum(N2) * log(lambda2) -
              sum(N1+N2) * log(lambda)
              )
    return(l)
  }

MEMLoglikeProfile <-  function(lambda1, Y)
  {
    op <- optimize(MEMLoglike, interval=c(0,2), maximum=T, lambda1=lambda1, Y=Y)
    return(op$objective)
  }

MEMIC <- function(Y, lambda1.est)
  {
    epsilon <- 0.01
    if(lambda1.est>0)
      {
        l.inf <- 1e-12
        l.sup <- lambda1.est
        l <- (l.inf+l.sup)/2
        maxloglike <- MEMLoglikeProfile(lambda1.est, Y)
        curloglike <- MEMLoglikeProfile(l, Y)
        while( (abs(maxloglike - curloglike -1.35) > epsilon )  & ( abs(l.inf-l.sup) > 1e-10) )
          {
                                        #print(c(l.inf, l.sup) )
            if( abs(maxloglike - curloglike) > (1.35 +epsilon) )
              {
                l.inf <- l
                l <- (l.inf+l.sup)/2
                curloglike <- MEMLoglikeProfile(l, Y)
              }else
            {
              l.sup <- l
              l <- (l.inf+l.sup)/2
              curloglike <- MEMLoglikeProfile(l, Y)
            }
          }
        l1 <- l
      }else
    {
      l1 <- 0
    }

    if(lambda1.est>0)
      {
        l.inf <- lambda1.est
      }else
    {
      l.inf <- 1e-12
      maxloglike <- MEMLoglikeProfile(l.inf, Y)
    }
    l.sup <- 1.8
    l <- (l.inf+l.sup)/2
    curloglike <- MEMLoglikeProfile(l, Y)
    while( (abs(maxloglike - curloglike -1.92) > epsilon ) & ( abs(l.inf-l.sup) > 1e-10) )
      {
        #print(c(l.inf, l.sup) )
        if( abs(maxloglike - curloglike) > (1.92 +epsilon) )
          {
            l.sup <- l
            l <- (l.inf+l.sup)/2
            curloglike <- MEMLoglikeProfile(l, Y)
          }else
        {
            l.inf <- l
            l <- (l.inf+l.sup)/2
            curloglike <- MEMLoglikeProfile(l, Y)
        }
      }
    l2 <- l
    return(c(l1,l2))
  }
