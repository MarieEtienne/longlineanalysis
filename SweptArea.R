computeSweptArea <- function( Y, ZoneCoverage=1)
  {
    return(mean( Y$N1/ (Y$P*Y$N)))
  }
