"brierScore.fitBMAgamma0" <-
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

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")
 
 L <- 1:nObs

 if (!is.null(dates) && !is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

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
 bsLogistic <- apply(MAT[L,,drop=FALSE]^2, 2, mean)

# BMA Brier Scores

  if (!any(is.na(WEIGHTS <- fit$weights))) {

    for (i in L) {
    
       f <- ensembleData[i,]
     
       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f
        
       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs,
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi)*fit$prob0coefs,
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
#            bma = 1 - safeDiv(bsBMA,bsClimatology))

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting,
            logistic = bsLogistic,
            bma = bsBMA)

}

