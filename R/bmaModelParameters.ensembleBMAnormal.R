"bmaModelParameters.ensembleBMAnormal" <-
function(fit, dates = NULL, ...) 
{
 if (is.null(dates)) dates <- names(fit$dateTable)

 I <- match(dates, names(fit$dateTable), nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for specified dates")

 list(weights = fit$weights[,I], 
      biasCoefs = fit$biasCoefs[,,I], 
      sd = if (is.null(dim(sd))) fit$sd[I] else fit$sd[,I],
      model = "normal")
}

