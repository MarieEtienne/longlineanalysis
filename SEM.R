###########################################################################
## Method 1  : Non linear regression as proposed by Somerton and Kikkawa 
###########################################################################
## script to perform the Non linear regression of the exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## Y the dataframe which stored the data, lomgline in this code
##############################################################################


## functions to be minimized to achieve non linear regression
## to stabilize optimisation procedure and to assure non negatvity of lambda
## working in log scale
ecart.exponential.model <- function(lambda.log, C, N, P)
  {
    lambda <- exp(lambda.log)
    return(sqrt(sum((C-N*(1-exp(-lambda*P)))^2)))
  }

functionToOptim <- function(lambda.log, C, N, P)
  {
        lambda <- exp(lambda.log)
	nb.esp <- ncol(C)
	contrib <- rep(0,nb.esp)	
	for (i in 1:nb.esp)
	{
		contrib[i] <- sum(
                   (C[,i]/N-lambda[i]/sum(lambda)*
                              (1-exp(-sum(lambda)*P)))^2
                            ) 
	}
	return( sqrt(sum(contrib) ))
  }



## main function for non linear regression estimations
AnalysisNLR <- function( subdata, tune1 = 3)
{ 
  ### Global estimation
  ## lambda.tot
  C.tot <- with(subdata, N-N0)
  nb.col <- dim(subdata)[2]
  Capt <- subdata[,6:(nb.col-1)]
  nb.species <- dim(Capt)[2]

  lambda.cur.log <- log(0.2)
  value.cur <- ecart.exponential.model(lambda.cur.log, C=C.tot, N=subdata$N, P=subdata$P)
  iter <- 0
  ok <- optimize(ecart.exponential.model, interval=c(-20,1), C=C.tot, N=subdata$N, P=subdata$P)
  
  
  
  lambda.cur.log <- log ( rep(exp(ok$minimum)/nb.species,nb.species) )
  value.cur <- functionToOptim(lambda.cur.log,  C=Capt, N=subdata$N, P=subdata$P)
  iter=0
  res <- nlm(functionToOptim, lambda.cur.log, C=Capt, N=subdata$N, P=subdata$P)
  #print(res$estimate)  
  return(list(lambda=exp(res$estimate), value=res$minimum, grad=res$gradient))  
}


## main function for non linear regression estimations 2 species only
AnalysisNLR2 <- function( subdata, tune1 = 3)
{

  if(!is.CLonglineData(subdata))
    {
      stop("Parameter ", subdata," of function AnalysisNLR2 is not of class CLonglineData")
    }
  else
    {
      
### Global estimation
      ## lambda.tot
      C.tot <- with(subdata, N-Nb)
      
      lambda.cur.log <- log(0.2)
      value.cur <- ecart.exponential.model(lambda.cur.log, C=C.tot, N=subdata$N, P=subdata$P)
      iter <- 0
      ok <- optimize(ecart.exponential.model, interval=c(-20,1), C=C.tot, N=subdata$N, P=subdata$P)
      
### Estimation species by species
      ## the function to use with package nleqslv
      ## for one value of lambda it returns the value of the system

      nb.col <- dim(subdata)[2]
      if(sum(subdata$Ne)>0)
        {
          Capt <- matrix(c(subdata$N1,subdata$N2, subdata$Ne), ncol=3)
        }else
      {
        Capt <- matrix(c(subdata$N1,subdata$N2), ncol=2)
      }
      nb.species <- dim(Capt)[2]
      
      lambda.cur.log <- log( rep(exp(ok$minimum) / nb.species,nb.species) )
      value.cur <- functionToOptim(lambda.cur.log,  C=Capt, N=subdata$N, P=subdata$P)

      res <- nlm(functionToOptim, lambda.cur.log, C=Capt, N=subdata$N, P=subdata$P)
                                        #print(res$estimate)  
	  if( sum(subdata$Ne)==0)
	  {
		res$estimate <- c(res$estimate,NA)
		res$minimum <- c(res$value,NA)
		res$gradient <- c(res$gradient,NA)
	  }
      return(list(lambda=exp(res$estimate), value=res$minimum, grad=res$gradient))  
    }
}


 
