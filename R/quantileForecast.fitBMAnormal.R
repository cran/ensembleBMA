`quantileForecast.fitBMAnormal` <-
function(fit, ensembleData, quantiles = 0.5, dates=NULL, ...) 
{
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")
 
 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS))) {

    SD <- if (!is.null(dim(fit$sd))) {
             fit$sd
          }
          else rep(fit$sd, nForecasts)

    for (i in 1:nObs) {
      
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       MEAN <- apply(rbind(1, f)*fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply(quantiles,quantBMAnormal,WEIGHTS=W,
                       MEAN=MEAN[!M],SD=SD[!M])
    }
 }

 Q
}

