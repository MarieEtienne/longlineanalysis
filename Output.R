## we is Empty if treated Empty hooks in the model
## NoEmpty if empty hooks are treated with non target species
ExploitSimu <- function(dir, nsimu, we ="Empty")
  {
    setwd(dir)

    file.init <- paste("SimuStudy", nsimu, ".init", sep="")
    file.out <-  paste("SimuStudy", nsimu, ".out", sep="")

    res <- read.table(file.out, skip=2, sep=",", header=T)
    theta <- ReadParameters(file.init)
    ##print(theta)
    
    ##selection of lines
    ind <- which(res$X==we)
    res <- res[ind,]
    res2 <- res

    if( theta$type=="without"  | we=="NoEmpty")
     {
       threshold <- 1
     }    else
      {
        threshold <- 0
      }

    M <- nrow(res)
    M2 <- M
    select <- sapply(1:M,  function(x) {  if( sum(is.na(res[x,])) > threshold ) {res2 <<- res2[-c(x),]; M2<<-M2-1} } )
    res <- res2
    rm(res2,M2)
    
    ## Estimation
    mle.l <- c(mean(res$mle.lambda1), mean(res$mle.lambda2))
    bayes.l <- c(mean(res$bayes.lambda1), mean(res$bayes.lambda2))
    nlr.l <- c(mean(res$nlr.lambda1), mean(res$nlr.lambda2))

    ##check if swept is in the simulations
    sw.bool= (sum(names(res)=="swept")>0)
    if(sw.bool) swept.l <- mean(res$swept)
    
    ## Bias Computation
    bias.mle.lY <- NA
    bias.bayes.lY <- NA
    bias.nlr.lY <- NA
    bias.swept.lY <- NA

    
    bias.mle.lY <- mean((res$mle.lambda1-theta$lambda1)/theta$lambda1)
    bias.bayes.lY <- mean((res$bayes.lambda1-theta$lambda1)/theta$lambda1)
    bias.nlr.lY <- mean((res$nlr.lambda1-theta$lambda1)/theta$lambda1)
    if(sw.bool) bias.swept.lY <- mean((res$swept-theta$lambda1)/theta$lambda1)

    ## CV Computation
    cv.mle.lY <- NA
    cv.bayes.lY <- NA
    cv.nlr.lY <- NA
    cv.swept.lY <- NA
    
    cv.mle.lY <- sd(res$mle.lambda1)/mean(res$mle.lambda1)
    cv.bayes.lY <- sd(res$bayes.lambda1)/mean(res$bayes.lambda1)
    cv.nlr.lY <- sd(res$nlr.lambda1)/mean(res$nlr.lambda1)
    if(sw.bool) cv.swept.lY <- sd(res$swept)/mean(res$swept)

    ## IC Computation
    ic.mle.lY <- NA
    ic.bayes.lY <- NA
    ic.nlr.lY <- NA
    
    ic.mle.lY <- sort(res$mle.lambda1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.bayes.lY <- sort(res$bayes.lambda1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.nlr.lY <- sort(res$nlr.lambda1)[c(floor(M*0.05), ceiling(M*0.95))]
    ic.swept <-  sort(res$swept)[c(floor(M*0.05), ceiling(M*0.95))]
    
##    par(mfcol=c(3,1))
##    r1 <- range(res$mle.lambda1, res$mbayes.lambda1, res$nlr.lambda1)
##    hist(res$mle.lambda1, xlim=r1)
##    hist(res$bayes.lambda1, xlim=r1)
##    hist(res$nlr.lambda1, xlim=r1)


    return(list(bias=c(bias.mle.lY, bias.bayes.lY, bias.nlr.lY, bias.swept.lY), cv=c(cv.mle.lY, cv.bayes.lY, cv.nlr.lY, cv.swept.lY),
             ic=matrix(c(ic.mle.lY, ic.bayes.lY, ic.nlr.lY, ic.swept), byrow=T, ncol=2), nlr.l=nlr.l, mle.l=mle.l, bayes.l=bayes.l )
           )
  }


DrawTable <- function(table.to.draw, col.lab, row.lab, scale.col,val.x , val.y, title.to.draw="", width=6, height=4,file.out="fig")
  {
    n.col <- ncol(table.to.draw) 
    n.row <- nrow(table.to.draw) 

    layout(matrix(c(1,2,0,3), ncol=2), width=c(8,2), height=c(1,8))
    par(mar=c(1,1,1,1))
    plot(x=-10,xlim=c(0,8), ylim=c(0,2), col.lab=0,axes=FALSE, col.main=1, col.axis=0)
    text(4,1,title.to.draw, col=1)
    
    
    par(mar=c(1,1,1,1))
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
    
    dev.copy2pdf(file=paste(file.out,".pdf",sep=""), height=height, width=width)
    dev.copy2eps(file=paste(file.out,".eps",sep=""), height=height, width=width)

   layout(matrix(c(1), ncol=1), width=c(8), height=c(8))
   }

