
#setwd("E:/MarieEtienne/Longlines") ## directory with the data
setwd("C:/Users/Marie/Rockfish/Longline")

source("RCode/LongLineJags.R")
source("RCode/LongLineMethod1.R")
source("RCode/LongLineExponentialModel.R")
source("RCode/SweptArea.R")
source("RCode/Utils.R")

##Prepare the data
## create the data for Bugs
## read the data from file longlineData.csv
longline <- read.table("DataPropreLongline_noRCA.csv", skip=1, header=T, sep=",")
longline <- longline[-298,] ## error in this line too much catch regarding the number of hooks.


head(longline)
#replace NA by 0
longline[is.na(longline)] <- 0

##compute baited hooks
longline <- data.frame(Year=longline$Year,
                       AREA=longline$DFO_STAT_AREA_CODE,
                       P=longline$soaktime,
                       N=longline$Nb_hooks_deployed,
                       N0=with(longline,Nb_hooks_deployed-Nb_yelloweye_caught- Nb_dogfish_caught - Nb_invertebrate_scavengers_caught -  Nb_allotherspecies_caught - Nb_inanimate_objects - Nb_hooks_empty),
                       N1=(longline$Nb_yelloweye_caught),
                       sp2=with(longline, Nb_dogfish_caught),
                       sp3=with(longline, Nb_invertebrate_scavengers_caught),
                       sp4=with(longline,  Nb_allotherspecies_caught),
                       sp5=with(longline,  Nb_inanimate_objects),
                       Nempty = with(longline, Nb_hooks_empty),
                       N2=with(longline, Nb_dogfish_caught + Nb_invertebrate_scavengers_caught + Nb_allotherspecies_caught + Nb_inanimate_objects)
      )


head(longline)

## Select Year and area
Area <- sort(unique(longline$AREA))
Year <- sort(unique(longline$Year))


col.title <-paste("area", "year", "empty/noempty", "mle$lambda1","mle$lambda2" , "bayes$lambda1", "bayes$lambda2", "bayes$Ic1_0.0025", "bayes$Ic1_0.0975",  "bayes$Ic2_0.0025","bayes$Ic2_0.0975" , "nlr$lambda1","nlr$lambda2","nlr$lambdaempty", "swept", sep=",")
Add2File(col.title,"ResEstimation_noRCA.txt", append=FALSE)

for(area in Area)
{
  for(year in Year)
  {
    print(paste("Area=", area,", Year=", year,sep=""))
    ind <- with(longline, which(AREA==area & Year==year))
    if(sum(longline[ind,]$N1)>0)
    {

      subdata <- longline[ind,]

      ##empty hooks are considered as a species
      nlr <- AnalysisNLR2(subdata)
      mle <- AnalysisMLE(subdata)
      bayes <- AnalysisBayes(subdata, length.MC.chain=10000)
      swept <- computeSweptArea(subdata)
      
      ##empty hooks are considered as other
      subdata <- longline[ind,]
      subdata$N <- subdata$N
      subdata$N2 <- subdata$N2+subdata$Nempty
      subdata$Nempty <- rep(0,nrow(subdata))

      nlr.empty <- AnalysisNLR2(subdata)
      mle.empty <- AnalysisMLE(subdata)
      bayes.empty <- AnalysisBayes(subdata, length.MC.chain=10000)

      
      all <- as.vector(unlist(c(area, year, "empty",  mle$lambda, bayes$lambda, bayes$Ic1, bayes$Ic2, nlr$lambda, swept)))
      Add2File(all, file="ResEstimation_noRCA.txt")
      all.empty <- as.vector(unlist(c(area, year, "noempty", mle.empty$lambda, bayes.empty$lambda, bayes.empty$Ic1, bayes.empty$Ic2, 
                        nlr.empty$lambda, swept )))
      Add2File(all.empty,file="ResEstimation_noRCA.txt")
    }
  } 
}
