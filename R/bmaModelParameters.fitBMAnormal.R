"bmaModelParameters.fitBMAnormal" <-
function(object, ...) 
{
 list(weights = object$weights, 
      biasCoefs = object$biasCoefs, 
      sd = object$sd,
      model = "normal")
}

