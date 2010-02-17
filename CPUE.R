################################################
## Function to compute CPUE for longline data
################################################


CPUE <- function( Y, ZoneCoverage=1)
  {
    if( ! is.CLonglineData(Y))
      {
        stop("First argument is not of class CLonglineData")
      }
    else
      {
        with(Y,
             {
               return(sum(N1) / sum(P*N) )
             }
           )
      }
  }
