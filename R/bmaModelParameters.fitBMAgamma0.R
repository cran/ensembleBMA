"bmaModelParameters.fitBMAgamma0" <-
function(object, ...) 
{
 list(weights = object$weights, 
      popCoefs = object$popCoefs,
      biasCoefs = object$biasCoefs, 
      varCoefs = object$varCoefs,
      transformation = object$transformation,
      inverseTransformation = object$inverseTransformation,
      model = "gamma0")
}

