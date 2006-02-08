"controlBMAgamma0" <-
function(maxIter = Inf, eps = sqrt(.Machine$double.eps), 
         nEsteps = 1, transformation = function(x) x^(1/3),
         inverseTransformation = function(x) x^3,
         start = list(varCoefs = NULL, weights = NULL)) 
{

 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, eps = eps, nEsteps = nEsteps,
      transformation = transformation, 
      inverseTransformation = inverseTransformation,
      start = start) 
}

