# longline
#' longline is used to create a longline object. 
#' 
#' @param fact1 a factor vector, for which which indices will be computed
#' @param fact2 an opional  factor vector which should be accounted for in the 
#' definition of the indices
#' @param nt a numeric vector of the number of individuals of target species
#'  caught
#' @param nnt a numeric vector of the number of individuals of non target 
#' species caught
#' @param nb a numeric vector of  the number of still baited hooks
#' @param ne an optional vector of  the number of empty hooks, if NULL, 
#' then empty hooks are considered as missing information
#' @param s a numeric vector of soaktime
#' 
#' @export
#' @return an object of class longline 
#
#' @examples
#' N <- 20
#' n <- rep(200,N)
#' s <- rep(100, N)
#' dataSim <- rmultinom(n=N, size=n, prob=c(.3, .1, .2, .4) )
#' testLongline <- longline(fact1 = as.factor(sample(2004:2006, size=N,
#' replace=TRUE)),  nb=dataSim[1,], nt=dataSim[,2], nnt=dataSim[,3], ne=dataSim[,4], s=s)


longline <- function(fact1, fact2=NULL,  nt, nnt, nb, ne=NULL, s)
{
  ne.missing <- F
  nLines <- length(fact1)
  if(is.null(ne)){
    ne.missing <- T
    ne=rep(NA,nLines)
  }
  if( nLines == length(nt) &
      nLines == length(nnt) &
      nLines == length(nb) & 
      nLines == length(ne) & 
      nLines == length(s) )
    {
      if(ne.missing)
        n=nt+nnt+nb
      else
      {
        ne[is.na(ne)] <- 0
        n=nt+nnt+nb+ne
      }
      Fact1   <- as.factor(fact1)
      if(!is.null(fact2))
        Fact2   <- as.factor(fact2)
      NFact1  <- nlevels(Fact1)
      if(!is.null(fact2))
        NFact2  <- nlevels(Fact2)
      NData  <- nLines
      nb[is.na(nb)] <- 0
      ## same soaktime  within fact1,fact2 ?
      cv.soaktime <- sd(s,na.rm = T)  /mean(s, na.rm = T)
      sames <- (cv.soaktime < 0.01)
       
      res <- list(Fact1=Fact1, NFact1=NFact1, NData=NData, 
                  Ne=ne, N=n, NT=nt, NNT=nnt, Nb=nb, S=s, sameS=sames)        
      if(!is.null(fact2))
        res <- c(res, list(NFact2=NFact2, Fact2=Fact2))
    }
  else
    {
      print("Lengths of the vectors differ")
      return (-1)
    }
  class(res) <- c("longline", "list")
  return(res)
}



# is.longline
#' is.longline is used to checked wheter or not the object is of class longline
#' 
#' @param x any object to te tested
#' @export
#' @return a boolean, TRUE if x is of class longline, False otherwise
#
#' @examples
#' N <- 20
#' n <- rep(200,N)
#' s <- rep(100, N)
#' dataSim <- rmultinom(n=N, size=n, prob=c(.3, .1, .2, .4) )
#' testLongline <- longline(fact1 = as.factor(sample(2004:2006, size=N,
#' replace=TRUE)), nb=dataSim[,1], nt=dataSim[,2], nnt=dataSim[,3], ne=dataSim[,4], s=s)
#' is.longline(testLongline)
#' is.longline(x=2)
is.longline <- function(x)
{
  return( inherits(x, "longline") )
}


# ExtractLongline
#' ExtractLongline is a longline object extracts from the argument, which matches the condition fact1 and fact2
#' 
#' @param x a object of class longline
#' @param fact1 the value of Fact1 to be used for extraction
#' @param fact2 the value of Fact2 to be used for extraction
#' @param ind if provided, only the individuals ind will be extracted whatever fact2 and  fact1 are

#' 
#' @export
#' @return a longline object with only observations matching Fact1=fact1 and/or  optionnaly Fact2=fact2
#
#' @examples
#' N <- 20
#' n <- rep(200,N)
#' s <- rep(100, N)
#' dataSim <- rmultinom(n=N, size=n, prob=c(.3, .1, .2, .4) )
#' testLongline <- longline(fact1 = as.factor(sample(2004:2006, size=N,
#' replace=TRUE)), nb=dataSim[1,], nt=dataSim[2,], nnt=dataSim[3,], ne=dataSim[4,], s=s)
#' ex <- ExtractLongline(testLongline, fact1="2004")
#' print(ex)
ExtractLongline <- function(x, fact1=NULL, fact2=NULL, ind=NULL)
{  
  if(is.null(ind))
  {
   if(is.null(fact1)&is.null(fact2))
    ind <- 1:(x$NData) 
   if( is.null(fact1) & (!is.null(fact2)))
     ind <- which(x$Fact2%in%fact2) 
   if( !is.null(fact1) & (is.null(fact2)))
     ind <- which(x$Fact1%in%fact1) 
   if( (!is.null(fact1)) & (!is.null(fact2)))
     ind <- which( (x$Fact1%in%fact1) & (x$Fact2%in%fact2)) 
  }
  
  res<- lapply(names(x), function(d){
     switch(d, 
            NData=length(ind),
            NFact1=nlevels(droplevels(x$Fact1[ind])),
            Fact1=droplevels(x$Fact1[ind]),
            NFact2=nlevels(droplevels(x$Fact2[ind])),
            Fact2=droplevels(x$Fact2[ind]),
            Ne=x$Ne[ind],
            N=x$N[ind],
            NT=x$NT[ind],
            NNT=x$NNT[ind],
            Nb=x$Nb[ind],
            S=x$S[ind],
            sameS = x$sameS
     )
   })
  names(res)=names(x)
  class(res) <- c("longline", "list")
  return(res)
}


# ComputeCPUE.longline
#' ComputeCPUE.longline computes classical CPUE according to Fact1 and Fact2
#' 
#' @param x a longline object
#' @export
#' @return classical CPUE  computed for every level of Fact1 and Fact2
#
#' @examples
#' data(longlineEx)
#' ll1 <- longline(fact1=longlineEx$Year, fact2=longlineEx$Area,  nt=longlineEx$NT, nnt=longlineEx$NNT, nb=longlineEx$Nb, ne=longlineEx$Ne, s=longlineEx$soaktime)
#' ComputeCPUE.longline(ll1)
#' ll2 <- longline(fact1=longlineEx$Year, nt=longlineEx$NT, nnt=longlineEx$NNT, nb=longlineEx$Nb, ne=longlineEx$Ne, s=longlineEx$soaktime)
#' ComputeCPUE.longline(ll2)
ComputeCPUE.longline<- function(x)
{
  
  if(with(x,  exists("Fact2"))){
    x$FComb <- paste0(x$Fact1, "-", x$Fact2 )
  }  else{
    x$FComb <- x$Fact1
  }
  x$FComb <- as.factor(x$FComb)
  res<- lapply(levels(x$FComb), function(d){
    ind <- which(x$FComb == d)
    ll <- ExtractLongline(x,  ind=ind)
    if(with(x,  exists("Fact2")))
    {
      ID  <- strsplit(d, split = "-")[[1]]
     prov <- data.frame(Fact1=ID[1],  Fact2=ID[2], CPUE=sum(ll$NT) / sum(ll$S*ll$N))
    }else   {
      prov <- data.frame(Fact1=d,   CPUE=sum(ll$NT) / sum(ll$S*ll$N))
    }
    return(prov)
    })
  res <- Reduce('rbind', res)
  
  return(res)
}
