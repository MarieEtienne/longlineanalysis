args <- commandArgs(TRUE)

area <- type.convert(args[1])
year <- type.convert(args[2])


#area=12
#year=2003

setwd("/home/metienne/Documents/Recherche/RockFish/Longlines")

file.in <- paste("SimuA", area,"Y",year,".txt",sep="")
res <- read.table(file.in)

hist.mle.par <- hist(res[,3], freq=F, main=paste("Area ", area, ", Year ", year, ": MLE bias study", sep=""), xlab="lambda1", ylab="Density")
hist.bayes.par <- hist(res[,5], freq=F, main=paste("Area ", area, ", Year ", year, ": posterior mean estimator bias study", sep=""), xlab="lambda1", ylab="Density")

xlimits=range(c(hist.mle.par$breaks, hist.bayes.par$breaks ))
ylimits=range(c(0,hist.mle.par$intensities, hist.bayes.par$intensities ))

##pdf(paste("BiasMLEArea", area,"Year",year,".pdf", sep=""), height=8, width=8)

jpeg(paste("BiasMLEArea", area,"Year",year,".jpg", sep=""), quality=90, height=480, width=480)

hist.mle.par <- hist(res[,3], freq=F, main=paste("Area ", area, ", Year ", year, ": MLE bias study", sep=""), xlab="lambda1", ylab="Density", xlim=xlimits, ylim=ylimits)
lines(rep( res[1,1],2), c(0, max(hist.mle.par$intensities)), col=2, lwd=3)
dev.off()


##pdf(paste("BiasBayesArea", area,"Year",year,".pdf", sep=""), height=8, width=8)
jpeg(paste("BiasBayesArea", area,"Year",year,".jpg", sep=""), quality=90, height=480, width=480)
hist.par <- hist(res[,5], freq=F, main=paste("Area ", area, ", Year ", year, ": posterior mean estimator bias study", sep=""), xlab="lambda1", ylab="Density", xlim=xlimits, ylim=ylimits)
lines(rep( res[1,1],2), c(0, max(hist.bayes.par$intensities)), col=2, lwd=3)
dev.off()
