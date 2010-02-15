setwd("E:/MarieEtienne/Longlines/Simulations")
source("DrawSamples.R")

### Toy example for one sample
lambda1 <- 0.0002
lambda2 <- 0.0001
P <- 50
N <- 500
data <- DrawOneSample( lambda1, lambda2, 120, 150 )
lambda.est1 <- with(data, N1 / (   P * (N1+N2) ) *log (N / N0))
lambda.est1
with(data, print(c(N,N0,N1,N2)))
print(paste("lambda1.est=",lambda.est1,", lambda1.sim=", lambda1,sep=""))
lambda2.est <- with(data, N2/N1)*lambda.est1
print(paste("lambda2.est=", lambda2.est, ", lambda2.sim=", lambda2,sep=""))
##---------------------> It works quite well

### Toy example for several samples
## with the same soaktime
lambda1 <- 0.0001
lambda2 <- 0.0005
nb.lines <- 30
#P <- sample(120:200, nb.lines, replace=TRUE)
P <- rep(120,nb.lines)
N <- sample(90:220, nb.lines, replace = TRUE)
data <- DrawSamples( lambda1, lambda2, P, N )
lambda1.est <- with(data, sum(N1) / (sum (N1+N2) *P[1] ) *log (sum(N) / sum(N0) ))
lambda2.est <- with(data, sum(N2) / sum (N1) ) * lambda1.est
print(paste("lambda1.est=",lambda1.est,", lambda1.sim=", lambda1,sep=""))
print(paste("lambda2.est=", lambda2.est, ", lambda2.sim=", lambda2,sep=""))

## initial values for minimization algorithm
p <- mean(data$P)
lambda1.init <- with(data,sum(N1)/sum(N1+N2) * log( sum(N) /sum(N0) )/p)

#optimize(EquationNormaleLambda1, lower=0.00001, upper=0.1, Y=data)
#nlm(EquationNormaleLambda1, p=lambda1.init, Y=data)
algo.min <-  nlminb(lambda1.init, EquationNormaleLambda1, lower = 0.0000001, upper = 0.01, Y=data)
lambda1.est <- algo.min$par
lambda2.est <- with(data, sum(N2) / sum (N1) ) * lambda1.est
print(paste("lambda1.est=",lambda1.est,", lambda1.sim=", lambda1,sep=""))
print(paste("lambda2.est=", lambda2.est, ", lambda2.sim=", lambda2,sep=""))


### Toy example for several samples
## with different soaktime
lambda1 <- 0.0001
lambda2 <- 0.0004
nb.lines <- 30
P <- sample(100:200, nb.lines, replace=TRUE)
#P <- rep(120,nb.lines)
N <- sample(90:220, nb.lines, replace = TRUE)
data <- DrawSamples( lambda1, lambda2, P, N )

## initial values for minimization algorithm
p <- mean(data$P)
algo.min.init <- with(data,sum(N1)/sum(N1+N2) * log( sum(N) /sum(N0) )/p)
value.cur <- EquationNormaleLambda1(algo.min.init, data)
algo.min <-  nlminb(algo.min.init, EquationNormaleLambda1, lower = 0.0000001, upper = 0.01, Y=data)
iter=0
tune1=3


while(iter < 10 & algo.min$objective< value.cur)
{
  iter <- iter + 1
  algo.min.init<- algo.min$par
  value.cur <- algo.min$objective
  algo.min <-  nlminb(algo.min.init, EquationNormaleLambda1, lower = algo.min.init - tune1 * algo.min.init / 10 , upper =algo.min.init +  tune1 * algo.min.init/10, Y=data)
  print( paste(" Objective=", algo.min$objective, ", Par=", algo.min$par, sep=""))
 
}

lambda1.est <- algo.min$par
lambda2.est <- with(data, sum(N2) / sum (N1) ) * lambda1.est
print(paste("lambda1.est=",lambda1.est,sep=""))
print(paste("lambda2.est=", lambda2.est,sep=""))


