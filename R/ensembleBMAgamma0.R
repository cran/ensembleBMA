ensembleBMAgamma0 <-
function(ensembleData, trainingDays, dates = NULL, 
         control = controlBMAgamma0(), exchangeable = NULL)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 if (missing(trainingDays)) stop("trainingDays must be specified")

 call <- match.call()

 warmStart <- FALSE

 if (missing(trainingDays)) stop("trainingDays must be specified")

 ensMemNames <- ensembleMembers(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable, ensembleGroups(ensembleData),
                                  nForecasts)

# remove instances missing all forecasts, obs or dates

 M <- !dataNA(ensembleData)
 if (!all(M)) ensembleData <- ensembleData[M,]
 
 nObs <- nrow(ensembleData)
 if (!nObs) stop("no data")

 Dates <- as.character(ensembleValidDates(ensembleData))
 DATES <- sort(unique(Dates))

 julianDATES <- ymdhTOjul(DATES)
 incr <- min(1,min(diff(julianDATES))) ## incr may be fractional for hours

 forecastHour <- ensembleFhour(ensembleData)
 lag <- ceiling( forecastHour / 24 )

## dates that can be modeled by the training data (ignoring gaps)

 dates <- getDates( DATES, julianDATES, dates, trainingDays, lag, incr)
 juliandates <- ymdhTOjul(dates)
 nDates <- length(dates)

 if (is.null(control$prior)) {
# accomodates saved mean as an additional parameter
   prob0coefs <- array( NA, c(3, nForecasts, nDates),
                        dimnames = list(NULL, ensMemNames, dates))
 }
 else {
   prob0coefs <- array( NA, c(4, nForecasts, nDates),
                        dimnames = list(NULL, ensMemNames, dates))
 }
 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, dates))
 varCoefs <- array( NA, c(2, nDates), dimnames = list(NULL, dates))
 weights <- array( NA, c(nForecasts, nDates),
                      dimnames = list(ensMemNames, dates))

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- loglikelihood <- rep(0, nDates)
 names(nIter) <- names(loglikelihood) <- dates

 obs <- dataVerifObs(ensembleData)
 
 K <- 1:nForecasts

 L <- length(juliandates)
 twin <- 1:trainingDays

 cat("\n")
 l <- 0
 for(i in seq(along = juliandates)) {

    I <- (juliandates[i]-lag*incr) >= julianDATES
    if (!any(I)) stop("insufficient training data")

    j <- which(I)[sum(I)]

    if (j != l) {

      D <- as.logical(match(Dates, DATES[j:1], nomatch=0))
      nonz <- sum(obs[D] != 0)
      if (is.null(control$prior) && nonz < control$rainobs) {
        cat("insufficient nonzero training obs for date", dates[i], "...\n")
        next
      }
      
      twin <- (j+1) - (1:trainingDays)

      if (is.null(control$prior)) {

# attempt to extend the training period

         while (TRUE) {
           D <- as.logical(match(Dates, DATES[twin], nomatch=0))
           if (!any(D)) stop("this should not happen")
           d <- ensembleValidDates(ensembleData[D,])
#     if (length(unique(d)) != trainingDays) stop("wrong # of training days")
           nonz <- sum(obs[D] != 0)
           if (nonz >= control$rainobs) break
           if (min(twin) == 1) break
           twin <- max(twin):(min(twin)-1)
        }

        if (nonz < control$rainobs) {
          cat("insufficient nonzero training obs for date", dates[i], "...\n")
          next
        }

      }

      cat("modeling for date", dates[i], "...")

      kNA <- apply(ensembleForecasts(ensembleData[D,]), 2, 
                   function(x) all(is.na(x)))

      if (any(kNA)) {
        if (!is.null(x <- exchangeable)) x <- exchangeable[-K[kNA]]
        fit <- fitBMAgamma0(ensembleData[D,-K[kNA]], control = control,
                            exchangeable = x)
      }
      else {

        fit <- fitBMAgamma0(ensembleData[D,], control = control,
                            exchangeable = exchangeable)
      }

      l <- j ## last model fit
      trainTable[i] <- length(unique(Dates[D]))
      nIter[i] <- fit$nIter 
      loglikelihood[i] <- fit$loglikelihood
      if (warmStart) control$start$weights <- weights[,i]
      cat("\n")
     }
   else {
     trainTable[i] <- -abs(trainTable[i-1])
     nIter[i] <- -abs(nIter[i-1])
     loglikelihood[i] <- loglikelihood[i-1]
   }

   prob0coefs[,K[!kNA],i] <- fit$prob0coefs
   biasCoefs[,K[!kNA],i] <- fit$biasCoefs
   varCoefs[,i] <- fit$varCoefs
   weights[K[!kNA],i] <- fit$weights

 }

 structure(list(training = list(days=trainingDays,lag=lag,table=trainTable), 
                prob0coefs = prob0coefs, biasCoefs = biasCoefs, 
                varCoefs = varCoefs, weights = weights, nIter = nIter,
                exchangeable = exchangeable, power = fit$power,
		call = match.call()),
                forecastHour = forecastHour, 
                initializationTime = ensembleItime(ensembleData),
                class = c("ensembleBMAgamma0","ensembleBMA"))
}

