###########################################################################
## Method 1  : Non linear regression as proposed by Somerton and Kikkawa 
###########################################################################
## script to perform the Non linear regression of the exponetial model
## If the data are store in a file longlineData.csv,
## We may just use function Analysis with 
## Y the dataframe which stored the data, lomgline in this code
##############################################################################


## functions to be minimized to achieve non linear regression
  ecart.exponential.model <- function(lambda, C, N, P)
  {
    return(sqrt(sum((C-N*(1-exp(-lambda*P)))^2)))
  }

functionToOptim <- function(lambda, C, N, P)
  {
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

  lambda.cur <- 0.2
  value.cur <- ecart.exponential.model(0.2, C=C.tot, N=subdata$N, P=subdata$P)
  iter <- 0
  ok <- optimize(ecart.exponential.model, interval=c(0.000001,0.2), C=C.tot, N=subdata$N, P=subdata$P)
  
  while(iter < 10 & value.cur > ok$objective)
  {
    iter <- iter + 1
    value.cur <- ok$objective
    lambda.cur <- ok$minimum
    ok <- optimize(ecart.exponential.model, interval=c(lambda.cur-tune1 * lambda.cur/10,lambda.cur+tune1 * lambda.cur/10),  C=C.tot, N=subdata$N, P=subdata$P)    
    #print(paste( " Objective=", ok$objective, ", Minimum=", ok$minimum, sep="")) 
  }
	
  ### Estimation species by species
  # the function to use with package nleqslv
  # for one value of lambda it returns the value of the system

  nb.col <- dim(subdata)[2]
  
  Capt <- subdata[,6:(nb.col-1)]
  nb.species <- dim(Capt)[2]
  
  lambda.cur <- rep(lambda.cur/nb.species,nb.species)
  value.cur <- functionToOptim(lambda.cur,  C=Capt, N=subdata$N, P=subdata$P)
  iter=0
  res <- nlm(functionToOptim, lambda.cur, C=Capt, N=subdata$N, P=subdata$P)
  #print(res$estimate)  
  return(list(lambda=res$estimate, value=res$minimum, grad=res$gradient))  
}


## main function for non linear regression estimations 2 species only
AnalysisNLR2 <- function( subdata, tune1 = 3)
{ 
  ### Global estimation
  ## lambda.tot
  C.tot <- with(subdata, N-N0)

  lambda.cur <- 0.2
  value.cur <- ecart.exponential.model(0.2, C=C.tot, N=subdata$N, P=subdata$P)
  iter <- 0
  ok <- optimize(ecart.exponential.model, interval=c(0.000001,0.2), C=C.tot, N=subdata$N, P=subdata$P)
  
  while(iter < 10 & value.cur > ok$objective)
  {
    iter <- iter + 1
    value.cur <- ok$objective
    lambda.cur <- ok$minimum
    ok <- optimize(ecart.exponential.model, interval=c(lambda.cur-tune1 * lambda.cur/10,lambda.cur+tune1 * lambda.cur/10),  C=C.tot, N=subdata$N, P=subdata$P)    
    #print(paste( " Objective=", ok$objective, ", Minimum=", ok$minimum, sep="")) 
  }
	
  ### Estimation species by species
  # the function to use with package nleqslv
  # for one value of lambda it returns the value of the system

  nb.col <- dim(subdata)[2]
  if(sum(subdata$Nempty)>0)
    {
      Capt <- matrix(c(subdata$N1,subdata$N2, subdata$Nempty), ncol=3)
    }else
  {
    Capt <- matrix(c(subdata$N1,subdata$N2), ncol=2)
  }
  nb.species <- dim(Capt)[2]
  
  lambda.cur <- rep(lambda.cur/nb.species,nb.species)
  value.cur <- functionToOptim(lambda.cur,  C=Capt, N=subdata$N, P=subdata$P)
  iter=0
  res <- nlm(functionToOptim, lambda.cur, C=Capt, N=subdata$N, P=subdata$P)
  #print(res$estimate)  
  return(list(lambda=res$estimate, value=res$minimum, grad=res$gradient))  
}

##Example
#AnalysisNLR(longlineArea12Year2003)


 
