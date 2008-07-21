`modelParameters.ensembleBMAgamma0` <-
function(fit, dates = NULL, ...) 
{

 dateTable <- dimnames(fit$weights)[[2]]

 if (is.null(dates)) dates <- dateTable

 dates <- sort(unique(as.character(dates)))

 if (length(dates) > length(dateTable)) 
   stop("parameters not available for some dates")

 I <- match( dates, dateTable, nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for some dates")

 list(weights = fit$weights[,I], 
      prob0coefs = fit$prob0coefs[,,I],
      biasCoefs = fit$biasCoefs[,,I], 
      varCoefs = fit$varCoefs[,I],
      transformation = fit$transformation,
      inverseTransformation = fit$inverseTransformation,
      model = "gamma0")
}

