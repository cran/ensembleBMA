`modelParameters.fitBMAnormal` <-
function(fit, ...) 
{
 list(weights = fit$weights, 
      biasCoefs = fit$biasCoefs, 
      sd = fit$sd,
      model = "normal")
}

