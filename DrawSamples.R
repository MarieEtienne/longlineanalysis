### Draw samples according to the complete exponential model
DrawSamples <- function(lambda1, lambda2, P, N, area=12, year=2003)
{
  I<- length(N)
  
  Y <- data.frame(AREA= NULL, Year=NULL, N0 = NULL, N1 = NULL, N2 = NULL, N = NULL, P = NULL)
  for(i in 1:I){
    T1 <- rexp(N[i], rate= lambda1 )
    T2 <- rexp(N[i], rate= lambda2 )
    T <- pmin(T1, T2)

    ind.empty <- which(T>P[i])
    ind.species1 <- which(T1<P[i] & T2>T1)
    ind.species2 <- which(T2 < P[i] & T1 > T2)
    N0 = sum(T > P[i])
    N1 = sum(T1 < P[i] & T2 > T1)
    N2 = sum(T2 < P[i] & T1 > T2)
    Yrow <- data.frame(AREA=area, Year=year, N0 = N0, N1 = N1, N2 = N2, N = N[i], P = P[i])
    Y<- rbind (Y, Yrow)
  }

  Y <- as.CLonglineData(Y$Year, Y$AREA, n=Y$N, n1= Y$N1, n2=Y$N2, nb=Y$N-Y$N1-Y$N2, ne=rep(0, I), p=Y$P)

 return(Y)

}

### Draw samples according to the complete exponential and allow escapement
DrawSampleswEscape <- function(par.simu)
{
  
  sample.data <- with(par.simu,DrawSamples(lambda1, lambda2, rep(P, L) , rep(N, L), area=05, year=2020) )

  if(par.simu$type=="without")
  {
    sample.data$Nempty <- rep(0, nrow(sample.data))
  }
  else if( par.simu$type=="uniform")
  {
    n1.escape <- rbinom( nrow(sample.data), sample.data$N1, prob=par.simu$percent)
    n2.escape <- rbinom( nrow(sample.data), sample.data$N2, prob=par.simu$percent)
   sample.data$Nempty <- n1.escape + n2.escape
    	
    sample.data$N1 <- sample.data$N1 -n1.escape
    sample.data$N2 <- sample.data$N2 - n2.escape
  }
  else if (par.simu$type=="preferential")
  {
    n1.escape <- rbinom( nrow(sample.data), sample.data$N1, prob=par.simu$percent)
    n2.escape <- rbinom( nrow(sample.data), sample.data$N2, prob=par.simu$preffacteur*par.simu$percent)
    sample.data$Nempty <- n1.escape + n2.escape
    	
    sample.data$N1 <- sample.data$N1 -n1.escape
    sample.data$N2 <- sample.data$N2 - n2.escape
   }
  
  sample.data <- with( sample.data,
                      as.CLonglineData(YEAR, AREA, n=N, n1= N1, n2=N2, nb=N-N1-N2-Nempty, ne=Nempty, p=P)
                      )
##  print(sample.data)
  return(sample.data)
}

### Draw samples according to the complete exponential and allow escapement
### input is real design
DrawSampleswEscape2 <- function(lambda1, lambda2, Y, type="witouht", percent=0.1, pref=1)
{
  
  sample.data <- DrawSamples(lambda1, lambda2, Y$P , Y$N, area=05, year=2020) 

  if(type=="without")
  {
    sample.data$Nempty <- rep(0, nrow(sample.data))
  }
  else if( type=="uniform")
  {
    n1.escape <- rbinom( nrow(sample.data), sample.data$N1, prob=percent)
    n2.escape <- rbinom( nrow(sample.data), sample.data$N2, prob=percent)
   sample.data$Nempty <- n1.escape + n2.escape
    	
    sample.data$N1 <- sample.data$N1 -n1.escape
    sample.data$N2 <- sample.data$N2 - n2.escape
  }
  else if (type=="preferential")
  {
    n1.escape <- rbinom( nrow(sample.data), sample.data$N1, prob=percent)
    n2.escape <- rbinom( nrow(sample.data), sample.data$N2, prob=preffacteur*percent)
    sample.data$Nempty <- n1.escape + n2.escape
    	
    sample.data$N1 <- sample.data$N1 -n1.escape
    sample.data$N2 <- sample.data$N2 - n2.escape
   }
  
  sample.data <- with( sample.data,
                      as.CLonglineData(YEAR, AREA, n=N, n1= N1, n2=N2, nb=N-N1-N2-Nempty, ne=Nempty, p=P)
                      )
##  print(sample.data)
  return(sample.data)
}

