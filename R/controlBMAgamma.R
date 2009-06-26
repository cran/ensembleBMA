`controlBMAgamma` <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), 
         power = 1, start = list(varCoefs = NULL, weights = NULL)) 
{

 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, power = power, start = start) 
}

