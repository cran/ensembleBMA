"brierSkillScores.ensembleBMAnormal" <-
function(object, ensembleData, thresholds, ...) 
{
 inverseLogit <- function(x) {
# logit function safeguared against underflow and overflow
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  1/(1+exp(-x))
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  if (x >= log(.Machine$double.eps)) {
                    x <- exp(x)
                    x/(1+x)
                  }
                  else exp(x)
                }
                else 0
             }
            }

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 nObs <- sum(object$dateTable)

 ensDates <- ensembleDates(ensembleData)
 Dates <- as.numeric(ensDates)
 DATES <- sort(unique(Dates))

 K <- sapply(names(object$dateTable), function(d,D) 
                        which(d == as.character(D))[1],
                            D = ensDates)
 dates <- sort(Dates[K])
 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 nForecasts <- ensembleSize(ensembleData) 

 L <- as.logical(match(Dates, 
   dates[dates >= DATES[object$trainingRule$length + object$trainingRule$lag]],
   nomatch = 0))

 y <- ensembleVerifObs(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)
 x <- apply( ensembleData, 1, mean)

 MAT <- matrix( NA, nrow = length(y), ncol = length(thresholds))
 
 MAT <-  outer(y, thresholds, "<=")

# wrong
# bsClimatology <- apply((apply(MAT,2,mean) - MAT)^2, 2, mean)

 bsClimatology <- apply(sweep(MAT[L,], MARGIN = 2, FUN = "-", 
                        STATS = apply(MAT[L,],2,mean))^2, 2, mean)

 bsVoting <- apply((t(apply(ensembleData[L, ], 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean),
                 thresholds = thresholds)) - MAT[L,])^2, 2, mean)

 offset <- 1 - object$trainingRule$lag - (1:object$trainingRule$length)

 for (j in J) {
# logistic Brier Scores

    if (any(j + offset < 1)) next

    TrainSet <- as.logical(match(Dates, DATES[j+offset], nomatch = 0))

    logisticFit <- sapply( thresholds, 
            function(thresh, x, y) 
             glm((y <= thresh) ~ x, family = binomial(logit))$coef,
             x = x[TrainSet], y = y[TrainSet])

    logisticFit[2,][is.na(logisticFit[2,])] <- 0

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    MAT[I,] <- apply(logisticFit, 2, function(coefs,x) 
                      sapply(coefs[1] + coefs[2]*x, inverseLogit),
                      x = x[I]) - outer(y[I], thresholds, "<=")
 }

 bsLogistic <- apply(MAT[L,]^2, 2, mean)

 k <- 0
 for (j in J) {
# BMA Brier Scores

    k <- k + 1

    if (any(is.na(object$weights[,k]))) next
     
    SD <- if (!is.null(dim(object$sd))) object$sd[,k] else rep(object$sd[k], nForecasts)

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * object$biasCoefs[,,k], 2, sum)

       WEIGHTS <- object$weights[,k]

       MAT[i,] <- sapply( thresholds, normalBMAcdf, offset = 0, 
              WEIGHTS = WEIGHTS, MEAN = MEAN, SD = SD) -
                    (y[i] <= thresholds)

    }

 }

# locations at which forecasts are made (depends on training length and lag)

 bsBMA <- apply(MAT[L,]^2, 2, mean)
 
 safeDiv <- function(x,y) {
              yzero <- !y
              nz <- sum(yzero)
              result <- rep(NA, length(y))
              if (!nz) result <- x/y else result[!yzero] <- x[!yzero]/y[!yzero]
              result
            }  

 data.frame(thresholds = thresholds,
            ensemble = 1 - safeDiv(bsVoting,bsClimatology), 
            logistic = 1 - safeDiv(bsLogistic,bsClimatology),  
            bma = 1 - safeDiv(bsBMA,bsClimatology))
}

