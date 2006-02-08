"brierScore.ensembleBMAgamma0" <-
function(fit, ensembleData, thresholds, dates = NULL, popData = NULL, ...) 
{
 
 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 inverseLogit <- function(x) {
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

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)

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

 if (is.null(ensDates <- ensembleDates(ensembleData))) {
   if (length(dateTable) > 1) stop("date ambiguity")
   Dates <- rep(1,nObs)
   dates <- DATES <- 1
   L <- 1:nObs
 }
 else {
   if (!is.null(dates)) {
     L <- as.logical(match(dates, as.character(ensDates), nomatch=0))
     if (all(!L) || !length(L)) 
       stop("specified dates incompatible with ensemble data")
   }
   Dates <- as.numeric(ensDates)
   DATES <- sort(unique(Dates))
   L <- as.logical(match(as.character(ensDates), names(dateTable), nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   dates <- sort(unique(Dates[L]))
 }

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 nForecasts <- ensembleSize(ensembleData) 

 if (is.null(y <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

 ensembleData <- ensembleForecasts(ensembleData)
 x <- sapply(apply( ensembleData, 1, mean), fit$transformation)
 
 MAT <-  outer(y, thresholds, "<=")

# wrong
# bsClimatology <- apply((apply(MAT,2,mean) - MAT)^2, 2, mean)

 bsClimatology <- apply(sweep(MAT[L,,drop=FALSE], MARGIN = 2, FUN = "-", 
                        STATS = apply(MAT[L,,drop=FALSE],2,mean))^2, 2, mean)

 bsVoting <- apply((t(apply(ensembleData[L, ], 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean),
                 thresholds = thresholds)) - MAT[L,])^2, 2, mean)

 offset <- 1 - fit$trainingRule$lag - (1:fit$trainingRule$length)

 MAT <- matrix( NA, nrow = nObs, ncol = length(thresholds))

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

 bsLogistic <- apply(MAT[L,,drop=FALSE]^2, 2, mean)

 l <- 0
 for (j in J) {
# BMA Brier Scores

    l <- l + 1
    k <- K[l]

    if (any(is.na(WEIGHTS <- fit$weights[,k])))  next

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       VAR <- fit$varCoefs[1,k] + fit$varCoefs[2,k]*f
        
       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,,k], 2, sum)

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs[,,k],
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi)*fit$prob0coefs[,,k],
                                    2,sum), inverseLogit)
       }

       MAT[i,] <- sapply( sapply(thresholds,fit$transformation), 
                          gamma0BMAcdf, 
          WEIGHTS=WEIGHTS, PROB0=PROB0, MEAN=MEAN, VAR=VAR) -
                                                (y[i] <= thresholds)

    }

 }

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
#

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting, 
            logistic = bsLogistic,
            bma = bsBMA)
}

