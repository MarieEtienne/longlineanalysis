
#setwd("E:/MarieEtienne/Longlines") ## directory with the data
#setwd("C:/Users/Marie/Rockfish/Longline")

source("SvnRCode/LongLineJags.R")
source("SvnRCode/LongLineMethod1.R")
source("SvnRCode/LongLineExponentialModel.R")
source("SvnRCode/SweptArea.R")
source("SvnRCode/Utils.R")

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

file.out="ResEstimation_noRCA.txt"
Add2File(x=c("Area", "Year" ,"","mle.lambda1", "mle.lambda2", "nlr.lambda1", "nlr.lambda2","nlr.lambdae","swept"), file=file.out, append=FALSE)
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
      sem1<- SEM.MLE(subdata)
      mem1 <-MEM.MLE(subdata)
      #bayes <- AnalysisBayes(subdata, length.MC.chain=10000)
      sw <- CPUE(subdata)
      sem2 <- SEM.MLE(subdata, SEM=2)
      mem2 <- MEM.MLE(subdata, MEM=2)
      #bayes.empty <- AnalysisBayes(subdata, length.MC.chain=10000)      
      res <- as.vector(c(area, year, mem1, mem2, sem1,sem2, sw)
      Add2File(all, file="ResEstimation_noRCA.txt")
    }
  } 
}
