"bmaModelParameters.fitBMAgamma0" <-
function(fit, ...) 
{
 list(weights = fit$weights, 
      popCoefs = fit$popCoefs,
      biasCoefs = fit$biasCoefs, 
      varCoefs = fit$varCoefs,
      transformation = fit$transformation,
      inverseTransformation = fit$inverseTransformation,
      model = "gamma0")
}

