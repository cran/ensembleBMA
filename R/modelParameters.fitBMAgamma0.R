`modelParameters.fitBMAgamma0` <-
function(fit, ...) 
{
 list(weights = fit$weights, 
      popCoefs = fit$popCoefs,
      biasCoefs = fit$biasCoefs, 
      varCoefs = fit$varCoefs,
      power = fit$power,
      model = "gamma0")
}

