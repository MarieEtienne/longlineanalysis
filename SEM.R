###########################################################################
## Method 1  : Non linear regression as proposed by Somerton and Kikkawa 
###########################################################################
## script to perform the Non linear regression of the exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## Y the dataframe which stored the data, lomgline in this code
##############################################################################

SEMLoglikePartial <- function(lambdasans1.l, lambda1.l,  Y) 
{
  lambda.l <- c(lambda1.l,lambdasans1.l)                               
  return( SEMLoglike(lambda.l,  Y) )
}


SEMLoglike <- function(lambda.l,  Y) 
{
  if(sum(Y$Ne)>0)
    {
      C <- matrix(c(Y$N1,Y$N2, Y$Ne), ncol=3)
    }else
  {
    C <- matrix(c(Y$N1,Y$N2), ncol=2)
  }
  N=Y$N
  P=Y$P
  nb.esp <- ncol(C)
  lambda <- exp(lambda.l)
  lambdaplus <- sum(lambda)
  contrib <- rep(0,nb.esp)
  n <- nb.esp * nrow(C)

  mu <- 0*C
  for (i in 1:nb.esp)
    {
      mu[,i] <-  C[,i] - N * lambda[i] / lambdaplus *
        ( 1 - exp(-lambdaplus * P) )
    }
  mu <- as.vector(mu)
  return( -n/2 * ( 1 + log(  sum(mu^2)/n ) + log(2*pi) ) )
}

SEMAIC <-  function(Y) 
{
  lambda.l <- log(rep(1e-4,3))
  res <- optim(par=lambda.l, SEMLoglike, method="BFGS", control=list(fnscale=-1), Y=Y)
  if(sum(Y$Ne)>0)
    {
      nbPar <- 4
    }
  else
    {
      nbPar <- 3
    }
  return(-2*res$value + 2 * nbPar)
}
  
  
SEMLoglikeProfile <-  function(l1, Y, init)
  {
    l1.l <- log(l1)
    init.l <- log(init)
    op <- optim(par=init.l, SEMLoglikePartial,method="BFGS", control=list(fnscale=-1), lambda1.l=l1.l, Y=Y)
    return(op$value)
  }

SEMIC <- function(Y, lambda.est)
  {
    lambda1.est <- lambda.est[1]
    if(lambda1.est>0)
      {
        l.inf <- 1e-12
        l.sup <- lambda1.est
        init <- lambda.est[-1]
        epsilon <- lambda1.est /10
        l <- (l.inf+l.sup)/2
        maxloglike <- SEMLoglikeProfile(lambda1.est, Y, init)
        curloglike <- SEMLoglikeProfile(l, Y, init)
        iter=0
        while( (abs(maxloglike - curloglike -1.35) > epsilon )& (abs(l.sup-l.inf)> 1e-10) )
          {
            iter <- iter+1
                                        #print(c(l.inf, l.sup) )
            if( abs(maxloglike - curloglike) > (1.35 +epsilon) )
              {
                l.inf <- l
                l <- (l.inf+l.sup)/2
                curloglike <- SEMLoglikeProfile(l, Y, init)
              }else
            {
              l.sup <- l
            l <- (l.inf+l.sup)/2
              curloglike <- SEMLoglikeProfile(l, Y, init)
            }
          }
        l1 <- l
      }else
    {l1 <- 0}
    
    if(lambda1.est>0)
      {
        l.inf <- lambda1.est
      }else
    {l.inf <- 1e-12}
    l.sup <- min(lambda1.est*1e5,1.8)
    l <- (l.inf+l.sup)/2
    maxloglike <- SEMLoglikeProfile(lambda1.est, Y, init)
    curloglike <- SEMLoglikeProfile(l, Y, init)
    iter =0
    while( (abs(maxloglike - curloglike -1.92) > epsilon ) &  (abs(l.sup-l.inf)> 1e-10))
      {
        iter <- iter+1
        #print(c(l.inf, l.sup) )
        if( abs(maxloglike - curloglike) > (1.92 +epsilon) )
          {
            l.sup <- l
            l <- (l.inf+l.sup)/2
            curloglike <- SEMLoglikeProfile(l, Y, init)
          }else
        {
            l.inf <- l
            l <- (l.inf+l.sup)/2
            curloglike <- SEMLoglikeProfile(l, Y, init)
        }
      }
    l2 <- l
    return(c(l1,l2))
  }



## main function for non linear regression estimations
## 
AnalysisNLR <- function(subdata, tune1 = 3)
{
  if(!is.CLonglineData(subdata))
    {
      stop("Parameter ", subdata," of function AnalysisNLR2 is not of class CLonglineData")
    }
  else
    {
      
      nb.col <- dim(subdata)[2]
      if(sum(subdata$Ne)>0)                                                      
        {
          Capt <- matrix(c(subdata$N1,subdata$N2, subdata$Ne), ncol=3)
        }else
      {
        Capt <- matrix(c(subdata$N1,subdata$N2), ncol=2)
      }
      nb.species <- dim(Capt)[2]
      
      lambda.l <- log(rep(1e-4, nb.species))
      res <- optim(par=lambda.l, SEMLoglike, method="BFGS", control=list(fnscale=-1), Y=subdata)
    
                                        #print(res$estimate)  
	  if( sum(subdata$Ne)==0)
	  {
		res$par <- c(res$par,NA)
	  }
      return(list(lambda=exp(res$par), value=res$value))  
    }
}


 
