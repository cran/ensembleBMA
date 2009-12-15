controlBMAgamma0 <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), power = (1/3),
         rainobs = 10, init = list(varCoefs = NULL, weights = NULL)) 
{

 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, power = power, rainobs = rainobs,
      init = init) 
}

