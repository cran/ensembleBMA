"bmaModelParameters.ensembleBMAgamma0" <-
function(fit, dates = NULL, ...) 
{
 if (is.null(dates)) dates <- names(fit$dateTable)

 I <- match(dates, names(fit$dateTable), nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for specified dates")

 list(weights = fit$weights[,I], 
      prob0coefs = fit$prob0coefs[,,I],
      biasCoefs = fit$biasCoefs[,,I], 
      varCoefs = fit$varCoefs[,I],
      transformation = fit$transformation,
      inverseTransformation = fit$inverseTransformation,
      model = "gamma0")
}

