"quantileForecastBMA.fitBMAnormal" <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, ...) 
{
 
 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 nObs <- ensembleNobs(ensembleData)

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")
 
 L <- 1:nObs

 if (!is.null(dates) && !is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 if (!any(is.na(WEIGHTS <- fit$weights))) {

    SD <- if (!is.null(dim(fit$sd))) {
             fit$sd
          }
          else rep(fit$sd, nForecasts)

    for (i in L) {
      
       f <- ensembleData[i,]

       MEAN <- apply(rbind(1, f)*fit$biasCoefs, 2, sum)

       Q[i,] <- sapply(quantiles,normalBMAquant,WEIGHTS=WEIGHTS,
                       MEAN=MEAN,SD=SD)
    }
 }

 Q
}

