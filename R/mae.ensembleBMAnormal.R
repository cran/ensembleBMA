`mae.ensembleBMAnormal` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 weps <- 1.e-4
 
 if (!is.null(seed)) set.seed(seed)

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

 sampleMedian <- sampleMean <- predictiveMean <- rep(NA,nObs)
 names(sampleMedian) <- ensembleObsLabels(ensembleData)

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
          else {
            rep(fit$sd[k], nForecasts)
          }

    VAR <- SD*SD

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

       predictiveMean[i] <- sum( W * MEAN[!M])

       if (!is.null(nSamples)) {

## MAE via simulation

         MEAN <- MEAN[!M]
         SD <- SD[!M]

         SAMPLES <- sample((1:nForecasts)[!M],size=nSamples,
                            replace=TRUE,prob=W) 
         tab <- table(SAMPLES)

         Z <- tab == 0

         tab <- tab[!Z]

         MEAN <- MEAN[!Z]
         SD <- SD[!Z]

       SAMPLES <- apply(cbind(seq(along = tab), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN, SD = SD)

       SAMPLES <- as.vector(unlist(SAMPLES))

       sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
    }
  }
 }

## L <- which(!is.na(crpsSim))
## l <- length(L)

## maeCli <- mean(abs(obs - median(obs)))

## maeEns <- mean(abs(obs - apply(ensembleData, 1, median)))

 maeCli <- mean(abs(obs - mean(obs)))
 maeEnsMean <- mean(abs(obs - apply(ensembleData, 1, mean, na.rm = TRUE)))
 maeEnsMedian <- mean(abs(obs - apply(ensembleData, 1, median, na.rm = TRUE)))
 mseEnsMean <- sum((obs - apply(ensembleData, 1, mean, na.rm = TRUE))^2)/length(obs)
 mseEnsMedian <- sum((obs - apply(ensembleData, 1, median, na.rm = TRUE))^2)/length(obs)

 if (is.null(nSamples)) {
     maeBMAmedian <- mean(abs(obs - Q)) 
     maeBMAmean <- mean(abs(obs - predictiveMean)) 
     mseBMAmedian <- sum((obs - Q)^2)/length(obs)
     mseBMAmean <- sum((obs - predictiveMean)^2)/length(obs)
 }
 else {
   maeBMAmedian <- mean(abs(obs - sampleMedian))
   maeBMAmean <- mean(abs(obs - sampleMean))
     mseBMAmedian <- sum((obs - sampleMedian)^2)/length(obs)
     mseBMAmean <- sum((obs - sampleMean)^2)/length(obs)
 }

##c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)

A <- array( 0, c(2,2,2), 
  dimnames = list(c("mean", "median"), c("ensemble", "BMA"), c("mae", "mse")))

A[, , 1] <- matrix(c(maeEnsMean, maeEnsMedian, maeBMAmean, maeBMAmedian), 2, 2)
A[, , 2] <- matrix(c(mseEnsMean, mseEnsMedian, mseBMAmean, mseBMAmedian), 2, 2)

c(ensemble = A[2,1,1], BMA = A[2,2,1])
c(ensemble = maeEnsMedian, BMA = maeBMAmedian)
}

