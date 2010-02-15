### Script to plot the results of several methods of estimation
### If results are stored in longline_survey03_08_expon_mod_Res.csv
### Just Use function 
### PlotArea(area, Y=yellow.results, file=paste("Area",area,sep=""), type="ps")

#setwd("E:/MarieEtienne/Longlines") ## directory where the results file are stored

estimations <- read.table("OveralResults.csv",
    skip=1, header=T, sep=",") ## data frame with the results

yellow.results <- data.frame(Year=estimations$Year, Area=estimations$Area, 
                             rel.dens = estimations$rel.dens.yellow, 
                             rel.dens.with.R = estimations$lambdaYellow,
                             mle = estimations$mle.lambda1,
                             SweptArea=estimations$SweptArea,
                             bayes = estimations$b.lambda1,
                             bayes.ci.inf = estimations$IcLambcda1inf,
                             bayes.ci.sup = estimations$IcLambcda1sup,
                             simu.ci.inf = estimations$IcSimuInf,
                             simu.ci.sup = estimations$IcSimuSup,
                             boot.ci.inf = estimations$IcBootInf,
                             boot.ci.sup = estimations$IcBootSup
                            ) # renamed the data frame


## function to call to obtain graphs
## area is the number of the tudied area
## Y the dataframe with the stored results
## Couleur to be used for the points in the graph
## file name of the file to save the output graph without extension
## type = ps or pdf type the output file
PlotArea <- function(area, Y, couleur=c(2,3,4,5,6), file="", type="pdf",pos="bottomright")
{
  erreur <- 0
 if( file!="")
 {
    if(type=="pdf")
	{ 
        filename <- paste(file,".pdf",sep="")
	  pdf(filename,width=8, height=8)
       }else
	{
        if(type=="ps")
        { 
          filename <- paste(file,".ps",sep="")
	    postscript(filename,width=12, height=12)
        }else
        {
           print("Available output file are pdf or ps")
           erreur = -1
	  }
      }
  }
  if (erreur==0)
  {
    ind <- which(Y$Area==area)
    subdata <- Y[ind,]
    par(new=FALSE)
    r1 <- with(subdata, runif(30, min=0, max=6e-04))
    with(subdata,plot(Year,rel.dens, pch=19, axes =FALSE, cex=1.5, col=couleur[1], ylim=c(0, 6e-004),bty = "n", main=paste("Area ", area,sep="")))
    with(subdata,points(Year,rel.dens.with.R, pch=19, cex=1.3, col=couleur[2]))
    with(subdata,points(Year,mle, pch=19, cex=1.1, col=couleur[3]))
    with(subdata,points(Year,bayes, pch=19, cex=1.1, col=couleur[4]))

    nb.year <- with(subdata,length(Year))
    for(i in 1:nb.year)
    {
	with(subdata,print(c(bayes.ci.inf[i],bayes.ci.sup[i])) )
      with(subdata, lines(c(Year[i],Year[i]), c(bayes.ci.inf[i],bayes.ci.sup[i]), col=couleur[4], lwd=2))     
    }
    with(subdata, axis(side=2, at = pretty(range(r1) ) ))
    with(subdata, axis(side=1, at = 2003:2008) )

    par(new=TRUE) # to plot the next graphic on the same layer
    with(subdata,plot(Year,SweptArea, pch=19, cex=1.2, col=couleur[5],  axes = FALSE, ylim=c(0, 902),bty = "n", xlab = "", ylab = ""))
    with(subdata, axis(side=4, at = pretty(range(c(0,902)))))

    legend( pos, legend=c("Rel Density", "Rel density with R", "MLE", "bayes","Swept Area"),
    col=couleur, pch=rep(19, 5))
    if(file!="")
    {
      dev.off()
    }
  }
}



###Example
#area=17
#PlotArea(area, Y=yellow.results, file=paste("Area",area,sep=""), type="ps")
#PlotArea(area, Y=yellow.results, file=paste("Area",area,sep=""), type="pdf")



PlotIC <- function(area, Y, couleur=c(2,3,4,5,6), file="", type="pdf",pos="bottomright")
{
  erreur <- 0
 if( file!="")
 {
    if(type=="pdf")
	{ 
        filename <- paste(file,".pdf",sep="")
	  pdf(filename,width=8, height=8)
       }else
	{
        if(type=="ps")
        { 
          filename <- paste(file,".ps",sep="")
	    postscript(filename,width=12, height=12)
        }else
        {
           print("Available output file are pdf or ps")
           erreur = -1
	  }
      }
  }
  if (erreur==0)
  {
    ind <- which(Y$Area==area)
    subdata <- Y[ind,]
    par(new=FALSE)
    r1 <- with(subdata, runif(30, min=0, max=6e-04))
    with(subdata,plot(Year,rel.dens, pch=19, axes =FALSE, cex=1.5, col=couleur[1], ylim=c(0, 6e-004),bty = "n", main=paste("Area ", area,sep="")))
    with(subdata,points(Year,rel.dens.with.R, pch=19, cex=1.3, col=couleur[2]))
    with(subdata,points(Year,mle, pch=19, cex=1.1, col=couleur[3]))
    with(subdata,points(Year,bayes, pch=19, cex=1.1, col=couleur[4]))

    nb.year <- with(subdata,length(Year))
    for(i in 1:nb.year)
    {
	with(subdata,print(c(bayes.ci.inf[i],bayes.ci.sup[i])) )
        with(subdata, lines(c(Year[i],Year[i]), c(bayes.ci.inf[i],bayes.ci.sup[i]), col=couleur[4], lwd=2))     
    }
    with(subdata, axis(side=2, at = pretty(range(r1) ) ))
    with(subdata, axis(side=1, at = 2003:2008) )

    par(new=TRUE) # to plot the next graphic on the same layer
    with(subdata,plot(Year,SweptArea, pch=19, cex=1.2, col=couleur[5],  axes = FALSE, ylim=c(0, 902),bty = "n", xlab = "", ylab = ""))
    with(subdata, axis(side=4, at = pretty(range(c(0,902)))))

    legend( pos, legend=c("Rel Density", "Rel density with R", "MLE", "bayes","Swept Area"),
    col=couleur, pch=rep(19, 5))
    if(file!="")
    {
      dev.off()
    }
  }
}



