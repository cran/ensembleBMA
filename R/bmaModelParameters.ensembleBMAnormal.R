"bmaModelParameters.ensembleBMAnormal" <-
function(object, dates = NULL, ...) 
{
 if (is.null(dates)) dates <- names(object$dateTable)

 I <- match(dates, names(object$dateTable), nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for specified dates")

 list(weights = object$weights[,I], 
      biasCoefs = object$biasCoefs[,,I], 
      sd = if (is.null(dim(sd))) object$sd[I] else object$sd[,I],
      model = "normal")
}

