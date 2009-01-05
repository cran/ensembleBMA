`controlBMAgamma` <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), 
         nEsteps = 1, power = 1,
         start = list(varCoefs = NULL, weights = NULL)) 
{

 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, nEsteps = nEsteps,
      power = power, start = start) 
}

