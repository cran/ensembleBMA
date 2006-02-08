"quantileForecastBMA.ensembleBMAgamma0" <-
function(object, ensembleData, quantiles = 0.5, popData = NULL, ...) 
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

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 ensDates <- ensembleDates(ensembleData)

 Dates <- as.numeric(ensDates)
 DATES <- sort(unique(Dates))

 L <- as.logical(match(as.character(ensDates),
                 names(object$dateTable), nomatch = 0))

 K <- sapply(names(object$dateTable), function(d,D) 
                        which(d == as.character(D))[1],
                            D = ensDates)
 dates <- sort(Dates[K])
 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")
 
 nObs <- ensembleNobs(ensembleData)
 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsNames(ensembleData),as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

 k <- 0
 for (j in J) {

    k <- k + 1

    if (any(is.na(WEIGHTS <- object$weights[,k]))) next

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))
    
    for (i in I) {
    
       f <- ensembleData[i,]

       VAR <- object$varCoefs[1,k] + object$varCoefs[2,k]*f

       fTrans <- sapply(f, object$transformation)

       if (is.null(popData)) {
         POP <- sapply(apply(rbind( 1, fTrans, f==0) * object$popCoefs[,,k],
                                  2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         POP <- sapply(apply(rbind( 1, fTrans, popi) * object$popCoefs[,,k],
                                2,sum), inverseLogit)
       }

       MEAN <- apply(rbind(1, fTrans)*object$biasCoefs[,,k], 2, sum)

       Q[i,] <- sapply(quantiles,gamma0BMAquant,WEIGHTS=WEIGHTS,
                       POP=POP,MEAN=MEAN,VAR=VAR)
    }
 }

  apply(Q[L,,drop = FALSE], 2, object$inverseTransformation)
}

