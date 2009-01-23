`MAE.default` <-
function(fit, ensembleData, dates=NULL, ...) 
{
 weps <- 1.e-4

 matchITandFH(fit,ensembleData)

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]
 
## match specified dates with dateTable in fit

 dateTable <- dimnames(fit$weights)[[2]]

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

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 obs <- ensembleVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData) 

 Q <- as.vector(quantileForecast( fit, ensembleData, dates = dates))

 ensembleData <- ensembleForecasts(ensembleData)

## maeCli <- mean(abs(obs - median(obs)))

## maeEns <- mean(abs(obs - apply(ensembleData, 1, median)))

 maeCli <- mean(abs(obs - mean(obs)))

 maeEnsMean <- mean(abs(obs - apply(ensembleData, 1, mean, na.rm = TRUE)))
 maeEnsMedian <- mean(abs(obs - apply(ensembleData, 1, median, na.rm = TRUE)))
 mseEnsMean <- sum((obs - apply(ensembleData, 1, mean, na.rm = TRUE))^2)/length(obs)
 mseEnsMedian <- sum((obs - apply(ensembleData, 1, median, na.rm = TRUE))^2)/length(obs)

     maeBMAmedian <- mean(abs(obs - Q)) 
     mseBMAmedian <- sum((obs - Q)^2)/length(obs)

##c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)

 c(ensemble = maeEnsMedian, BMA = maeBMAmedian)
}

