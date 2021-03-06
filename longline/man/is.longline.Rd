% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longlineData.R
\name{is.longline}
\alias{is.longline}
\title{is.longline is used to checked wheter or not the object is of class longline}
\usage{
is.longline(x)
}
\arguments{
\item{x}{any object to te tested}
}
\value{
a boolean, TRUE if x is of class longline, False otherwise
}
\description{
is.longline is used to checked wheter or not the object is of class longline
}
\examples{
N <- 20
n <- rep(200,N)
s <- rep(100, N)
dataSim <- rmultinom(n=N, size=n, prob=c(.3, .1, .2, .4) )
testLongline <- longline(fact1 = as.factor(sample(2004:2006, size=N,
replace=TRUE)), nb=dataSim[,1], nt=dataSim[,2], nnt=dataSim[,3], ne=dataSim[,4], s=s)
is.longline(testLongline)
is.longline(x=2)
}

