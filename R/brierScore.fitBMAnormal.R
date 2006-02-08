"brierScore.fitBMAnormal" <-
function(fit, ensembleData, thresholds, dates = NULL, ...) 
{

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 nObs <- ensembleNobs(ensembleData)

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")
 
 L <- 1:nObs

 if (!is.null(dates) && is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

 nForecasts <- ensembleSize(ensembleData) 

 if (is.null(y <- ensembleVerifObs(ensembleData))) 
   stop("verification observations required")

 ensembleData <- ensembleForecasts(ensembleData)
 x <- apply( ensembleData, 1, mean)
 
 MAT <-  outer(y, thresholds, "<=")

 bsClimatology <- apply(sweep(MAT[L,,drop=FALSE], MARGIN = 2, FUN = "-", 
                        STATS = apply(MAT[L,,drop=FALSE],2,mean))^2, 2, mean)

 bsVoting <- apply((t(apply(ensembleData[L, ], 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean),
                 thresholds = thresholds)) - MAT[L,])^2, 2, mean)

 offset <- 1 - fit$trainingRule$lag - (1:fit$trainingRule$length)

 MAT <- matrix( NA, nrow = nObs, ncol = length(thresholds))

 bsLogistic <- apply(MAT[L,,drop=FALSE]^2, 2, mean)

# BMA Brier Scores

 if (!any(is.na(WEIGHTS <- fit$weights))) {
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd 
          }
         else rep(fit$sd, nForecasts)

    for (i in L) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)

       MAT[i,] <- sapply( thresholds, normalBMAcdf,
                         WEIGHTS = WEIGHTS, MEAN = MEAN, SD = SD) -
                         (y[i] <= thresholds)

    }

 }

# locations at which forecasts are made (depends on training length and lag)

 bsBMA <- apply(MAT[L,,drop=FALSE]^2, 2, mean)
 
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
            logistic = bsLogistic,
            bma = bsBMA)
}

