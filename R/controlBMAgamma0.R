"controlBMAgamma0" <-
function(maxIter = .Machine$integer.max, eps = sqrt(.Machine$double.eps), 
         nEsteps = 1, transformation = function(x) x^(1/3),
         inverseTransformation = function(x) x^3,
         start = list(varCoefs = NULL, weights = NULL)) 
{

 list(maxIter = maxIter, eps = eps, nEsteps = nEsteps,
      transformation = transformation, 
      inverseTransformation = inverseTransformation,
      start = start) 
}

