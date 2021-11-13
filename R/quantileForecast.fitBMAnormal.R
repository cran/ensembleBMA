quantileForecast.fitBMAnormal <-
function (fit, ensembleData, quantiles = 0.5, dates = NULL, ...) 
{
    weps <- 1e-04
    if (!is.null(dates)) 
        warning("dates ignored")
    M <- matchEnsembleMembers(fit, ensembleData)
    nForecasts <- ensembleSize(ensembleData)
    if (is.null(exchangeable <- ensembleGroups(ensembleData)))
      exchangeable <- 1:nForecasts  
    if (!all(M == 1:nForecasts)) 
        ensembleData <- ensembleData[, M]
    ensembleData <- ensembleData[!dataNA(ensembleData, observations = FALSE, 
        dates = FALSE), ]
    nObs <- dataNobs(ensembleData)
    Q <- matrix(NA, nObs, length(quantiles))
    dimnames(Q) <- list(dataObsLabels(ensembleData), as.character(quantiles))
    nForecasts <- ensembleSize(ensembleData)
    ensembleData <- ensembleForecasts(ensembleData)
    WEIGHTS <- fit$weights
    if (!all(Wmiss <- is.na(WEIGHTS))) {
        SD <- as.vector(fit$sd)
        if (fit$control$equalVariance) {
	  SD <- rep(SD,nForecasts)
	}
        if (length(SD) != nForecasts) 
            stop("length(SD) != nForecasts")
        for (i in 1:nObs) {
            f <- ensembleData[i, ]
            M <- is.na(f) | Wmiss
            MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)
            W <- WEIGHTS
            if (any(M)) {
                W <- W + weps
                W <- W[!M]/sum(W[!M])
            }
            Q[i, ] <- sapply(quantiles, quantBMAnormal, WEIGHTS = W, 
                MEAN = MEAN[!M], SD = SD[!M])
        }
    }
    Q
}
