################################################
## Function to compute CPUE for longline data
################################################


computeCPUE <- function( Y, ZoneCoverage=1)
  {
    if( ! is.CLonglineData(Y))
      {
        stop("First argument is not of class CLonglineData")
      }
    else
      {
        return(mean(ZoneCoverage * Y$N1/ (Y$P*Y$N)))
      }
  }
