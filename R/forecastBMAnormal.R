"forecastBMAnormal" <-
function(date, ensembleData, trainingRule = list(length = 30, lag = 2), 
         control = controlBMAnormal(), minCRPS = FALSE, quantiles = 0.5)
{

 if (!inherits(ensembleData, "ensembleData")) stop("not an ensembleData object")
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

 D <- as.logical(match(as.numeric(ensembleData$dates), DATES[I], nomatch = 0))

 fit <- fitBMAnormal(ensembleData[D,], control = control)

 if (minCRPS) {
      CRPSobjective <- function(SD) {
                crpsNormal(SD, fit$weights, fit$biasCoefs, ensembleData[D,])
            }
      if (control$equalVariance) {
        fit$sd <- optimize(CRPSobjective, interval = c(0, 6*fit$sd))$minimum
      }
      else {
        opt <- optim(fit$sd, CRPSobjective)
        if (!opt$convergence) {
          fit$sd <- opt$par
        }
        else {
          warning("CRPS could not be minimized")
        }
      }
 }

 SD <- if (!is.null(dim(fit$sd))) fit$sd else rep(fit$sd, nForecasts)

 D <- date == as.numeric(ensDates)

 I <- which(D)

 Q <- matrix(0, nrow = length(I), ncol = length(quantiles))
 if (length(quantiles) == 1 && quantiles == 0.5) {
   dimnames(Q) <- list(ensembleObsNames(ensembleData)[I], NULL)
  }
 else {
dimnames(Q) <- list(ensembleObsNames(ensembleData)[I], as.character(quantiles))
  }

 ensembleData <- ensembleForecasts(ensembleData)

 k <- 0
 for (i in I) {

    f <- ensembleData[i,]

    MEAN <- apply(rbind(1, f)*fit$biasCoefs, 2, sum)

    k <- k + 1

    Q[k,] <- sapply(quantiles,normalBMAquant,WEIGHTS=fit$weights,
                       MEAN=MEAN,SD=SD)
 }

 Q
}

