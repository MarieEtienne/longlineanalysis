setwd("/home/metienne/Documents/Recherche/RockFish/Longlines")


area=12
year=2007

sim <- read.table(paste("SimuA",area,"Y",year,".txt",sep=""), sep=" ", header=F, col.names=c("lambda1.true", "lambda2.true", "lambda1.mle", "lambda2.mle", "lambda1.bayes", "lambda2.bayes"))

av.bias.mle <- with(sim, mean((lambda1.mle -lambda1.true)/lambda1.true) )
print(av.bias.mle)

av.bias.bayes <- with(sim, mean((lambda1.bayes -lambda1.true)/lambda1.true) )
print(av.bias.bayes)
