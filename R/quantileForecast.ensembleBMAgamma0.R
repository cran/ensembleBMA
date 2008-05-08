`quantileForecast.ensembleBMAgamma0` <-
function(fit, ensembleData, quantiles=0.5, dates=NULL, ...) 
{
 weps <- 1.e-4

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 nObs <- nrow(ensembleData)

 if (!is.null(dates)) {
   K <- match(dates, names(fit$dateTable), nomatch=0)
   if (any(!K) || !length(K)) 
     stop("parameters not available for a specified date")
   dateTable <- fit$dateTable[K]
 }
 else {
   dateTable <- fit$dateTable
   K <- 1:length(dateTable)
  }

 ## match specified dates with dateTable in fit

 dateTable <- names(fit$nIter)

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable)) 
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K)) 
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

## match dates in data with dateTable
 if (is.null(ensDates <- ensembleDates(ensembleData))) {
   if (length(dates) > 1) stop("date ambiguity")
   Dates <- rep( dates, nObs)
 }
 else {
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }
 
 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {
    
    l <- l + 1
    k <- K[l]

    WEIGHTS <- fit$weights[,k]
    if (all(Wmiss <- is.na(WEIGHTS))) next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))
    
    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- fit$varCoefs[1,k] + fit$varCoefs[2,k]*f

       fTrans <- sapply(f, fit$transformation)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f==0)*fit$prob0coefs[,,k],
                                    2,sum), inverseLogit)

       MEAN <- apply(rbind(1, fTrans)*fit$biasCoefs[,,k], 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply(quantiles, quantBMAgamma0, WEIGHTS=W,
                       PROB0=PROB0[!M], MEAN=MEAN[!M], VAR=VAR[!M])
    }
 }

  apply(Q, 2, fit$inverseTransformation)
}

