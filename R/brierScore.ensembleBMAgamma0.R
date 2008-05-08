`brierScore.ensembleBMAgamma0` <-
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

 dateTable <- names(fit$nIter)

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable)) 
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K)) 
     stop("parameters not available for a some dates")

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
   y <- y[L]
 }

 nObs <- length(y)

 nForecasts <- ensembleSize(ensembleData) 

 ensembleData <- ensembleForecasts(ensembleData)
 x <- sapply(apply( ensembleData, 1, mean, na.rm = TRUE), fit$transformation)
 
 MAT <-  outer(y, thresholds, "<=")

 bsClimatology <- apply(sweep(MAT, MARGIN = 2, FUN = "-", 
                        STATS = apply(MAT,2,mean))^2, 2, mean)

 bsVoting <- apply((t(apply(ensembleData, 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean, na.rm = TRUE),
                 thresholds = thresholds)) - MAT)^2, 2, mean, na.rm = TRUE)

# avoid training data and apply to all data

logisticFit <- sapply( thresholds, 
            function(thresh, x, y) 
             glm((y <= thresh) ~ x,family=binomial(logit))$coef,
             x = x, y = y)

 logisticFit[2,][is.na(logisticFit[2,])] <- 0

 MAT <- apply(logisticFit, 2, function(coefs,x) 
                      sapply(coefs[1] + coefs[2]*x, inverseLogit),
                      x = x) - outer(y, thresholds, "<=")

 bsLogistic <- apply(MAT^2, 2, mean)

 MAT <- matrix( NA, nrow = nObs, ncol = length(thresholds))
 dimnames(MAT) <- list(NULL, as.character(thresholds))

 l <- 0

 for (d in dates) {
# BMA Brier Scores

    l <- l + 1
    k <- K[l]

    WEIGHTS <- fit$weights[,k]
    if (any(Wmiss <- is.na(WEIGHTS)))  next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       VAR <- fit$varCoefs[1,k] + fit$varCoefs[2,k]*f
        
       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,,k], 2, sum)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs[,,k],
                              2,sum), inverseLogit)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       MAT[i,] <- sapply( sapply(thresholds,fit$transformation), 
                          cdfBMAgamma0, 
          WEIGHTS=W, PROB0=PROB0[!M], MEAN=MEAN[!M], VAR=VAR[!M]) -
                                                (y[i] <= thresholds)

    }

 }

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
#

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting, 
            logistic = bsLogistic,
            bma = bsBMA)
}

