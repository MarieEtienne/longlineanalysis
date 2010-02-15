args <- commandArgs(TRUE)

Nsimu <- type.convert(args[1])
area <- type.convert(args[2])
year <- type.convert(args[3])
  

#setwd("E:/MarieEtienne/Longlines") ## directory with the data
setwd("/home/metienne/Documents/Recherche/RockFish/Longlines")


source("RCode/LongLineJags.R")
source("RCode/LongLineExponentialModel.R")
source("RCode/DrawSamples.R")
source("RCode/SweptArea.R")

## To compute Swept Area, define Zone Covergae
ZoneCoverage <- matrix(c(12, 13, 14, 15, 16, 17, 18,
                         797.08, 454.37, 144.31, 242.64, 262.55,  282.82, 217.18), ncol=2)



##Prepare the data
## create the data for Bugs
## read the data from file longlineData.csv
longline <- read.table("DataPropreLongline_noRCA.csv", skip=1, header=T, sep=",")
longline <- longline[-298,]
head(longline)
#replace NA by 0
longline[is.na(longline)] <- 0

##compute baited hooks
longline <- data.frame(longline,
    N0=with(longline, Nb_hooks_deployed-Nb_yelloweye_caught - Nb_dogfish_caught- Nb_invertebrate_scavengers_caught- Nb_allotherspecies_caught- Nb_inanimate_objects- Nb_hooks_empty),
    N1=with(longline,Nb_yelloweye_caught),
    N2=with(longline, Nb_dogfish_caught + Nb_invertebrate_scavengers_caught + Nb_allotherspecies_caught + Nb_inanimate_objects + Nb_hooks_empty),
    N= with(longline,Nb_hooks_deployed),
    P= with(longline,soaktime)
      )
                       
head(longline)

ind <- with(longline, which(AREA==area & Year==year))
subdata <- longline[ind,]  
nrow(subdata)
sum(subdata$N)/ nrow(subdata)
sum(subdata$N)


Y=subdata

SimulationAnalysis <- function(Y, M)
{
  t1 <- proc.time()
  mle <- AnalysisMLE(Y)
  estimates <- matrix(0, ncol=2, nrow=M)	
  bayes.estimates <- matrix(0, ncol=2, nrow=M)	
  swept <- numeric(M)

  file.out <- paste("SimuA", Y$AREA[1],"Y",Y$Year[1],".txt",sep="")
  lost <- 0
  for(i in 1:M)
    {
      data.sim <- DrawSamples(mle$lambda[1], mle$lambda[2], Y$P, Y$N, year=Y$Year[1], area=Y$AREA[1])
      while( sum(data.sim$N1)== 0 & lost < 1000)
        {
          data.sim <- DrawSamples(mle$lambda[1], mle$lambda[2], Y$P, Y$N)
          lost <- lost+1
        }
      bayes.sim <- AnalysisBayes(data.sim)
      mle.sim <- AnalysisMLE(data.sim)
      bayes.estimates[i,] <- c(bayes.sim$lambda1, bayes.sim$lambda2)
      estimates[i,] <- mle.sim$lambda
      swept[i] <- computeSweptArea(data.sim, ZoneCoverage[ZoneCoverage[,1]==data.sim$AREA[1],2])
      write(c(mle$lambda[1], mle$lambda[2],estimates[i,], bayes.estimates[i,]),ncolumns=6,file.out, append=TRUE)
    }
  ic=sort(estimates[,1])[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.bayes=sort(bayes.estimates[,1])[c(floor(M * 0.025), ceiling(M * 0.975))]
  ic.swept=sort(swept)[c(floor(M * 0.025), ceiling(M * 0.975))]

  print(proc.time()-t1)
  
  return(list(true.parameters=mle$lambda,
              pt.estimates=estimates,
              IC.lambda1=ic,
              IC.bayes=ic.bayes,
              IC.swept=ic.swept,
              lost=lost)
         )
}


print(paste("Area=", area,", Year=", year,sep=""))
ind <- with(longline, which(AREA==area & Year==year))
if(length(ind)>=1)
  {
    subdata <- longline[ind,]
    sim.study <- SimulationAnalysis(subdata, Nsimu)
    write(c(area, year,sim.study$true.parameters, sim.study$IC.lambda1,sim.study$IC.bayes,sim.study$IC.swept), ncolumns= 10, "SimStudy.txt", append=T, sep=",")
  }else
{
  print(paste("Dataset for year ", year," in area ", area,  "is  empty. Choose another dataset"))
}
