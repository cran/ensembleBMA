"bmaModelParameters.ensembleBMAgamma0" <-
function(object, dates = NULL, ...) 
{
 if (is.null(dates)) dates <- names(object$dateTable)

 I <- match(dates, names(object$dateTable), nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for specified dates")

 list(weights = object$weights[,I], 
      popCoefs = object$popCoefs[,,I],
      biasCoefs = object$biasCoefs[,,I], 
      varCoefs = object$varCoefs[,I],
      transformation = object$transformation,
      inverseTransformation = object$inverseTransformation,
      model = "gamma0")
}

