#######################################
##
## Definition of object ClonglineData
##
#######################################


as.CLonglineData <- function(year, area, n, n1, n2, nb, ne, p)
{
  n.lines <- length(year)
  if(n.lines == length(area) &
     n.lines == length(n) &
     n.lines == length(n1) &
     n.lines == length(nb) & 
     n.lines == length(ne) & 
     n.lines == length(p) )
    {
       res <- data.frame(YEAR=year, AREA=area, N=n, N1=n1, N2=n2, Nb=nb, Ne=ne, P=p)
    }
  else
    {
      print("Lengths of the vectors differ")
      return (-1)
    }
  class(res) <- c("CLonglineData", "data.frame")
  return(res)
}



is.CLonglineData <- function(x)
{
  return( inherits(x, "CLonglineData") )
}


