
## function to minimize in order to obtain the MLestimator for the exponential model, 
## with several longlines, with different soaktimes
EquationNormaleLambda1 <- function(lambda2, Y)
{
  cn <- 1 + with(Y, sum(N1) / sum(N2) )
  res <- with ( Y, - sum (P * N0)  + 
                           sum(  (N- N0) * P * exp( - lambda2  * cn *P ) / 
                                 ( 1 - exp( - lambda2  * cn *P ) ) 
                               )
              )
    return( log(res^2+1) )
}
