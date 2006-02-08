"quantileForecastBMA.ensembleBMAnormal" <-
function(object, ensembleData, quantiles = 0.5, ...) 
{

 ensDates <- ensembleDates(ensembleData)

 Dates <- as.numeric(ensDates)
 DATES <- sort(unique(Dates))

 L <- as.logical(match(as.character(ensDates),
                 names(object$dateTable), nomatch = 0))


 K <- sapply(names(object$dateTable), function(d,D) 
                        which(d == as.character(D))[1],
                            D = ensDates)
 dates <- sort(Dates[K])
 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 nObs <- ensembleNobs(ensembleData)
 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsNames(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 k <- 0
 for (j in J) {
    k <- k + 1

    if (any(is.na(WEIGHTS <- object$weights[,k]))) next

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))


    SD <- if (!is.null(dim(object$sd))) object$sd[,k] else rep(object$sd[k], nForecasts)

    for (i in I) {
    
       f <- ensembleData[i,]

       MEAN <- apply(rbind(1, f)*object$biasCoefs[,,k], 2, sum)

       Q[i,] <- sapply(quantiles,normalBMAquant,
                       WEIGHTS=WEIGHTS,MEAN=MEAN,SD=SD)
    }
 }

 Q[L,, drop = FALSE]
}

