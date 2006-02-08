"quantileForecastBMA.fitBMAnormal" <-
function(object, ensembleData, quantiles = 0.5, ...) 
{

 nObs <- ensembleNobs(ensembleData)
 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsNames(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 SD <- if (!is.null(dim(object$sd))) {
          object$sd
       }
       else rep(object$sd, nForecasts)

 WEIGHTS <- object$weights

 for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       MEAN <- apply(rbind(1, f)*object$biasCoefs, 2, sum)

       Q[i,] <- sapply(quantiles,normalBMAquant,WEIGHTS=WEIGHTS,
                       MEAN=MEAN,SD=SD)
 }

 Q
}

