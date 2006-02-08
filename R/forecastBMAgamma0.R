"forecastBMAgamma0" <-
function(date, ensembleData, trainingRule = list(length = 30, lag = 2), 
         control = controlBMAgamma0(), quantiles = 0.5, popData = NULL)
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
   if (dim(popData) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 call <- match.call()

 nObs <- ensembleNobs(ensembleData)

 ensDates <- ensembleDates(ensembleData)

 if (is.character(date)) {
   i <- which(date == as.character(ensDates))[1]
   date <- as.numeric(ensDates)[i]
 }

 D <- as.logical(match(date, as.numeric(ensDates), nomatch = 0))
 if (!D) stop("no match for date in data")

 nForecasts <- ensembleSize(ensembleData)
 
 uDates <- unique(as.numeric(ensDates))
 DATES <- sort(uDates)
 nDates <- length(DATES)
 seqDates <- seq(from = min(DATES), to = max(DATES), by = 1)
 if (length(DATES) != length(seqDates)) warning("some dates are missing")

 trainingRule <- as.list(trainingRule)
 startLoc <-  date - (trainingRule$length + trainingRule$lag) + 1
 if (startLoc < 1) 
   stop("training length and lag inconsistent with forecast date")

 I <- startLoc + (1:trainingRule$length) - 1

 D <- as.logical(match(as.numeric(ensDates), DATES[I], nomatch = 0))

 fit <- fitBMAgamma0(ensembleData[D,], control = control)

 D <- date == as.numeric(ensDates)

 I <- which(D)

 Q <- matrix(0, nrow = length(I), ncol = length(quantiles))
 if (length(quantiles) == 1 || quantiles == 0.5) {
   dimnames(Q) <- list(ensembleNobs(ensembleData)[I], NULL)
  }
 else {
   dimnames(Q) <- list(ensembleNobs(ensembleData)[I], as.character(quantiles))
  }

 ensembleData <- ensembleForecasts(ensembleData)

 k <- 0
 for (i in I) {

    f <- ensembleData[i,]

    VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f

    fTrans <- sapply(f, fit$transformation)

    MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

    if (is.null(popData)) {
      POP <- sapply(apply(rbind( 1, fTrans, f==0) * fit$popCoefs,
                                2,sum), inverseLogit)
    }
    else {
      popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
      POP <- sapply(apply(rbind( 1, fTrans, popi) * fit$popCoefs,
                                2,sum), inverseLogit)
    }

    k <- k + 1
  
    Q[k,] <- sapply(quantiles,gamma0BMAquant,WEIGHTS=fit$weights,
                         POP=POP, MEAN=MEAN, VAR=VAR)
 }

 Q
}

