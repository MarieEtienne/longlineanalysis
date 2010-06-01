#setwd("C:\\Users\\Marie\\Rockfish\\Longline\\WorkingDirectory1")
file.in <- "SimuStudy2.init"
file.res <- "SimuStudy2.out"

simu <- read.table(file.res, skip=3, header=F, sep=",")
ind.empty <- which( simu[,1] == "Empty" )
ind.non.empty <- which( simu[,1] == "NoEmpty" )

par(mfcol=c(2,1))
plot(simu[ind.empty,2],simu[ind.non.empty,2])
abline(0,1)
plot(simu[ind.empty,4],simu[ind.non.empty,4])
abline(0,1)

par.simu <- ReadParameters(file.in)
n.simu <- par.simu$n.simu
bias.mle <- ( simu[,2] - par.simu$lambda1 ) / par.simu$lambda1
bias.bayes <- ( simu[,4] - par.simu$lambda1 ) / par.simu$lambda1
par(mfcol=c(2,2))
hist( bias.mle[ind.empty] )
hist( bias.bayes[ind.empty] )
hist( bias.mle[ind.non.empty] )
hist( bias.bayes[ind.non.empty] )

   mean( bias.mle[ind.empty] )
 mean( bias.mle[ind.non.empty] )
 mean( bias.bayes[ind.empty] )
 mean( bias.bayes[ind.non.empty] )

sd( simu[ind.empty,2] )/mean( simu[ind.empty,2] )