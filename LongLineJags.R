##############################################################
## Method 3  : Bayesian estimation 
##############################################################
## script to perform the Bayesian estimation of the complete exponetial model
## using BRugs package
##############################################################################
require(coda)
require("rjags")

AnalysisBayes <- function(subdata, init=c(0.0001, 0.001), thin=50, burnin=2000, length.MC.chain=5000)
{
  if(!is.CLonglineData(subdata))
    {
      stop("Parameter 1 of function AnalysisBayes is not of class CLonglineData")
    }
  else
    {
      if(sum(subdata$N1 == 0 & subdata$N2==0 )>0)
        {
          return(list(lambda=c(NA,NA),
                      Ic1=c(NA,NA),
                      Ic2=c(NA,NA)
                      )
                 )
          
        }
      else
        {
        
          area <- subdata$AREA[1]
          year <- subdata$YEAR[1]
          N <-  subdata$N
          P=subdata$P
          L=nrow(subdata) 
          
          ##init values
          
          jags.fit <- jags.model(file= "LonglineJagsModel.txt",
                                 data=list("N0"=subdata$Nb, "N1"=subdata$N1,
                                   "N2"=subdata$N2,
                                   "Nidentified"=subdata$N1+subdata$N2,
                                   "N"=N, "P"=P, "L"=L),
                                 inits=list("lambda[1]"=init[1], "lambda[2]"=init[2]),
                                 n.chains=1)
          
          mcmc.data <- coda.samples(jags.fit,c("lambda[1]", "lambda[2]"), n.iter=2000)
          stat <- summary(mcmc.data)
          
          return(list(lambda=c(stat$statistic[1,1],stat$statistic[2,1]),
                      Ic1=stat$quantiles[1,c(1,5)],
                      Ic2=stat$quantiles[1,c(1,5)]
                      )
                 )
        }
    }
}
