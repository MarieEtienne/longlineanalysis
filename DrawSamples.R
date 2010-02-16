### Draw samples according to the complete exponential model
DrawSamples <- function(lambda1, lambda2, P, N, p1=0, p2=0, area=12, year=2003)
{
  I<- length(N)
  lambda    = lambda1 + lambda2
  Nb        = rbinom(I, N, prob=exp(- lambda * P) )
  N1PlusNE1 = rbinom(I, N-Nb, prob=lambda1/lambda)
  N2PlusNE2 = N-Nb-N1PlusNE1
  N1        = rbinom(I, N1PlusNE1, prob=(1-p1))
  NE1       = N1PlusNE1 - N1
  N2        = rbinom(I, N2PlusNE2, prob=(1-p2))
  NE2       = N2PlusNE2 - N2
  NE        = NE1 + NE2
  Y <- as.CLonglineData(year=rep(year,I), area=rep(area, I), n=N, n1= N1, n2=N2, nb=Nb, ne=NE, p=P)
 return(Y)
}

### Draw samples according to the complete exponential and allow escapement
DrawSampleswEscape <- function(par.simu)
{
  if(par.simu$type=="without")
  {
    sample.data <- with(par.simu,DrawSamples(lambda1, lambda2, rep(P, L) , rep(N, L), area=05, year=2020) )
  }
  else if( par.simu$type=="uniform")
  {
    sample.data <- with(par.simu,DrawSamples(lambda1, lambda2, rep(P, L) , rep(N, L), p1=par.simu$percent, p2=par.simu$percent,  area=05, year=2020) )
  }
  else if (par.simu$type=="preferential")
  {
    sample.data <- with(par.simu,DrawSamples(lambda1, lambda2, rep(P, L) , rep(N, L), p1=par.simu$percent, p2=par.simu$preffacteur*par.simu$percent,  area=05, year=2020) )
  }
  return(sample.data)
}

### Draw samples according to the complete exponential and allow escapement
### input is real design
DrawSampleswEscape2 <- function(lambda1, lambda2, Y, type="without", percent=0.1, pref=1)
{
  if(type=="without")
  {
   sample.data <- DrawSamples(lambda1, lambda2, Y$P , Y$N, area=05, year=2020) 
  }
  else if( type=="uniform")
  {
    sample.data <- DrawSamples(lambda1, lambda2, Y$P , Y$N, p1=percent, p2=percent, area=05, year=2020) 
  }
  else if (type=="preferential")
  {
    sample.data <- DrawSamples(lambda1, lambda2, Y$P , Y$N, p1=percent, p2=pref*percent, area=05, year=2020) 
   }
  return(sample.data)
}

