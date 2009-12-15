controlBMAnormal <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), equalVariance = TRUE, 
         biasCorrection = c("regression", "additive", "none"),
         init = list(sd = NULL, weights = NULL))  
{
 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, equalVariance = equalVariance,
      biasCorrection = biasCorrection[1], init = init)
}

