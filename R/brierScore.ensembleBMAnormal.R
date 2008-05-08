`brierScore.ensembleBMAnormal` <-
function(fit, ensembleData, thresholds, dates = NULL, ...) 
{
 weps <- 1.e-4

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts, obs, or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 if (is.null(y <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

#nObs <- length(y)
 nObs <- ensembleNobs(ensembleData)

## match specified dates with dateTable in fit

 dateTable <- names(fit$nIter)

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable)) 
     stop("parameters not available for some dates")

   K <- match(  dates, dateTable, nomatch=0)

   if (any(!K) || !length(K)) 
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

 }

 dateTable

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
   obs <- obs[L]
   nObs <- length(obs)
 }

 ensembleData <- ensembleForecasts(ensembleData)
 
 MAT <-  outer(y, thresholds, "<=")

 bsClimatology <- apply(sweep(MAT, MARGIN = 2, FUN = "-", 
           STATS = apply(MAT,2,mean))^2, 2, mean, na.rm = TRUE)

 bsVoting <- apply((t(apply(ensembleData, 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean, na.rm = TRUE),
                 thresholds = thresholds)) - MAT)^2, 2, mean, na.rm = TRUE)

 MAT <- matrix( NA, nObs, length(thresholds))
 dimnames(MAT) <- list(NULL, as.character(thresholds))

 l <- 0
 for (d in dates) {
# BMA Brier Scores

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

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,k], 2, sum)

       W <- WEIGHTS
       if (any(M)) { 
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       MAT[i,] <- sapply( thresholds, cdfBMAnormal,
                         WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M]) -
                         (y[i] <= thresholds)

    }

 }

# locations at which forecasts are made (depends on training length and lag)

 bsBMA <- apply(MAT^2, 2, mean, na.rm = TRUE)
 
 safeDiv <- function(x,y) {
              yzero <- !y
              nz <- sum(yzero)
              result <- rep(NA, length(y))
              if (!nz) result <- x/y else result[!yzero] <- x[!yzero]/y[!yzero]
              result
            }  

# data.frame(thresholds = thresholds,
#            ensemble = 1 - safeDiv(bsVoting,bsClimatology), 
#            logistic = 1 - safeDiv(bsLogistic,bsClimatology),  
#            bma = 1 - safeDiv(bsBMA,bsClimatology))

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting,
            bma = bsBMA)
}

