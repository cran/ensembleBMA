`cdf.ensembleBMAnormal` <-
function(fit, ensembleData, values, dates = NULL, ...) 
{
 weps <- 1.e-4

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

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
   L <- as.logical(match(Dates, dates, nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData), as.character(values)) 

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0 
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    WEIGHTS <- fit$weights[,k]
     
    if (all(Wmiss <- is.na(WEIGHTS))) next
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd[,k] 
          }
          else rep(fit$sd[k], nForecasts)

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,k], 2, sum)

       M <- is.na(f) | Wmiss

       W <- WEIGHTS

       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       CDF[i,] <- sapply( values, cdfBMAnormal,
                          WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M]) 

    }

 }

 CDF
}

