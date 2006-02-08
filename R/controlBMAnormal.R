"controlBMAnormal" <-
function(maxIter = .Machine$integer.max, eps = sqrt(.Machine$double.eps),
         equalVariance = TRUE, 
         biasCorrection = c("regression", "additive", "none"),
         start = list(sd = NULL, weights = NULL)) 
{
 list(maxIter = maxIter, eps = eps, equalVariance = equalVariance,
      biasCorrection = biasCorrection[1], start = start)
}

