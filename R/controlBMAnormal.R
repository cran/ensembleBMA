"controlBMAnormal" <-
function(maxIter = Inf, eps = sqrt(.Machine$double.eps), equalVariance = TRUE, 
         biasCorrection = c("regression", "additive", "none"),
         start = list(sd = NULL, weights = NULL))  
{
 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, eps = eps, equalVariance = equalVariance,
      biasCorrection = biasCorrection[1], start = start)
}

