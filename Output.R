ExploitSimu <- function(dir, nsimu)
  {
    file.init <- paste(dir,"/SimuStudy", nsimu, ".init", sep="")
    file.out <-  paste(dir,"/SimuStudy", nsimu, ".out", sep="")

    print(file.init)
    res <- read.table(file.out, skip=2, sep=",", header=T)
    theta <- ReadParameters(file.init)
    ##print(theta)
    
    ## removes bad results with NA values in the estimates
    M <- nrow(res)
    listena=rep(T,M)
    select <- sapply(M:1,  function(x) {  if( sum(is.na(res[x,])) > 2 ) {  listena[x]=F } } )
    res = res[listena,] 
    ## Estimation
    mem1 <- c(mean(res$mem1.l1), mean(res$mem1.l2))
    mem2 <- c(mean(res$mem2.l1), mean(res$mem2.l2))
    sem1 <- c(mean(res$sem1.l1), mean(res$sem1.l2))
    sem2 <- c(mean(res$sem2.l1), mean(res$sem2.l2))
    sem1.g <- c(mean(res$sem1.g.l1), mean(res$sem1.g.l2))
    sem2.g <- c(mean(res$sem2.g.l1), mean(res$sem2.g.l2))
    cpue <- mean(res$sw)

    ##check if swept is in the simulations
    sw.bool= (sum(names(res)=="swept")>0)
    if(sw.bool) swept.l <- mean(res$swept)
    
    ## Bias Computation
    bias.mem1 <- NA
    bias.mem2 <- NA
    bias.sem1 <- NA
    bias.sem2 <- NA
    bias.sem1.g <- NA
    bias.sem2.g <- NA
    bias.cpue <- NA

   
    bias.mem1 <- mean((res$mem1.l1-theta$lambda1)/theta$lambda1)
    bias.mem2 <- mean((res$mem2.l1-theta$lambda1)/theta$lambda1)
    bias.sem1 <- mean((res$sem1.l1-theta$lambda1)/theta$lambda1)
    bias.sem2 <- mean((res$sem2.l1-theta$lambda1)/theta$lambda1)
    bias.sem1.g <- mean((res$sem1.g.l1-theta$lambda1)/theta$lambda1)
    bias.sem2.g <- mean((res$sem2.g.l1-theta$lambda1)/theta$lambda1)
    bias.cpue <- mean((res$sw-theta$lambda1)/theta$lambda1)
   
    ## CV Computation
    cv.mem1 <- NA
    cv.mem2 <- NA
    cv.sem1 <- NA
    cv.sem2 <- NA
    cv.sem1.g <- NA
    cv.sem2.g <- NA
    cv.cpue <- NA
    
    cv.mem1 <- sd(res$mem1.l1)/mean(res$mem1.l1)
    cv.mem2 <- sd(res$mem2.l1)/mean(res$mem2.l1)
    cv.sem1 <- sd(res$sem1.l1)/mean(res$sem1.l1)
    cv.sem2 <- sd(res$sem2.l1)/mean(res$sem2.l1)
    cv.sem1.g <- sd(res$sem1.g.l1)/mean(res$sem1.l1)
    cv.sem2.g <- sd(res$sem2.g.l1)/mean(res$sem2.l1)
    cv.cpue <- sd(res$sw)/mean(res$sw)
   
   
   
    ic.mem1 <- sort(res$mem1.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.mem2 <- sort(res$mem2.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.sem1 <- sort(res$sem1.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.sem2 <- sort(res$sem2.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.sem1.g <- sort(res$sem1.g.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.sem2.g <- sort(res$sem2.g.l1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.cpue <- sort(res$sw)[c(floor(M*0.05), ceiling(M*0.95))]
    
##    par(mfcol=c(3,1))
##    r1 <- range(res$mle.lambda1, res$mbayes.lambda1, res$nlr.lambda1)
##    hist(res$mle.lambda1, xlim=r1)
##    hist(res$bayes.lambda1, xlim=r1)
##    hist(res$nlr.lambda1, xlim=r1)

    return(list(
      bias=c(bias.mem1, bias.mem2, bias.sem1, bias.sem2, bias.sem1.g, bias.sem2.g, bias.cpue),
      cv=c(cv.mem1, cv.mem2, cv.sem1, cv.sem2, cv.sem1.g, cv.sem2.g, cv.cpue),
      ic=matrix(c(ic.mem1, ic.mem2, ic.sem1, ic.sem2, ic.sem1.g, ic.sem2.g, ic.cpue), byrow=T, ncol=2),
      value=res,
      theta= theta
      ))
  }


DrawTable <- function(table.to.draw, col.lab, row.lab, scale.col,val.x , val.y, title.to.draw="", width=6, height=4,file.out="fig", ptsize=12)
{
  n.col <- ncol(table.to.draw) 
  n.row <- nrow(table.to.draw) 

  pdf(file=paste(file.out,".pdf",sep=""), height=height, width=width, pointsize=ptsize)
    layout(matrix(c(1,2,0,3), ncol=2), width=c(8,2), height=c(1,8))
    par(mar=c(1,1,1,1))
    plot(x=-10,xlim=c(0,8), ylim=c(0,2), col.lab=0,axes=FALSE, col.main=1, col.axis=0)
    text(4,1,title.to.draw, col=1)
    ## Define color
    col <- matrix(sapply(abs(table.to.draw), function(x) { sum(x > scale.col)+1 }), ncol=4)
    couleur <-gray( 1- seq(0,1,1/(length(scale.col)+3) ) )
    ## Draw Layout 
    plot(x=-10,xlim=c(0,n.row+2),ylim=c(0,n.col+2),axes=FALSE, col.main=1,  col.lab=0)
    text(0.5,n.row+0.5,row.lab,col=1)
    lines(c(1,1),c(0,n.col+2),col=1)
    lines(c(0,0),c(0,n.col+2), col=1)
    lines(c(n.row+1,n.row+1),c(0,n.col+2), col=1)
    lines(c(0,n.row+1),c(n.col+2,n.col+2), col=1)
    lines(c(1,n.row+1),c(n.col+1,n.col+1), col=1)
    text(n.col/2+0.5, n.row+1.5, col.lab)
    lines(c(0,n.row+1),c(0,0), col=1)
    for( i in 1:n.row)
      {                                                                                                                            
        lines(c(0,n.col+1),c(i,i), col=1)
        text(0.5, n.row-i+0.5, val.x[i]) 
        for(j in 1:n.col)
          {
            if(i==1)  text(j+0.5, n.col+0.5, val.y[j])
            rect(j, n.row-i ,j+1, n.row-i+1, col=couleur[col[i,j]])
            text(j+0.5, n.row-i+0.5, round(table.to.draw[i,j]*100,1), col=1)
          }
      }
    ## legend
    par(mar=c(1,1,1,1))
    plot(x=-10,xlim=c(0,2),ylim=c(0,8),axes=FALSE, col.main=0, col.lab=0)
    step.height=7/length(scale.col)
    for( i in 1:length(scale.col))
      {
        rect(0, i*step.height, 1, (i+1)*step.height, col=couleur[i] )
        text(1.5,(i+0.5)* step.height, paste(scale.col[i]*100,"%",sep=" ") )
        
      }
  dev.off()

#    postscript(file=paste(file.out,".pdf",sep=""), height=height, width=width)
#    layout(matrix(c(1,2,0,3), ncol=2), width=c(8,2), height=c(1,8))
#    par(mar=c(1,1,1,1))
#    plot(x=-10,xlim=c(0,8), ylim=c(0,2), col.lab=0,axes=FALSE, col.main=1, col.axis=0)
#    text(4,1,title.to.draw, col=1)
#    ## Define color
#    col <- matrix(sapply(abs(table.to.draw), function(x) { sum(x > scale.col)+1 }), ncol=4)
#    couleur <-gray( 1- seq(0,1,1/(length(scale.col)+3) ) )
#    ## Draw Layout 
#    plot(x=-10,xlim=c(0,n.row+2),ylim=c(0,n.col+2),axes=FALSE, col.main=1,  col.lab=0)
#    text(0.5,n.row+0.5,row.lab,col=1)
#    lines(c(1,1),c(0,n.col+2),col=1)
#    lines(c(0,0),c(0,n.col+2), col=1)
#    lines(c(n.row+1,n.row+1),c(0,n.col+2), col=1)
#    lines(c(0,n.row+1),c(n.col+2,n.col+2), col=1)
#    lines(c(1,n.row+1),c(n.col+1,n.col+1), col=1)
#    text(n.col/2+0.5, n.row+1.5, col.lab)
#    lines(c(0,n.row+1),c(0,0), col=1)
#    for( i in 1:n.row)
#      {                                                                                                                            
#        lines(c(0,n.col+1),c(i,i), col=1)
#        text(0.5, n.row-i+0.5, val.x[i]) 
#        for(j in 1:n.col)
#          {
#            if(i==1)  text(j+0.5, n.col+0.5, val.y[j])
#            rect(j, n.row-i ,j+1, n.row-i+1, col=couleur[col[i,j]])
#            text(j+0.5, n.row-i+0.5, round(table.to.draw[i,j]*100,1), col=1)
#          }
#      }
#    ## legend
#    par(mar=c(1,1,1,1))
#    plot(x=-10,xlim=c(0,2),ylim=c(0,8),axes=FALSE, col.main=0, col.lab=0)
#    step.height=7/length(scale.col)
#    for( i in 1:length(scale.col))
#      {
#        rect(0, i*step.height, 1, (i+1)*step.height, col=couleur[i] )
#        text(1.5,(i+0.5)* step.height, paste(scale.col[i]*100,"%",sep=" ") )
#        
#      }
#      dev.off()
#
#    

   layout(matrix(c(1), ncol=1), width=c(8), height=c(8))
   }

