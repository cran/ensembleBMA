`ensembleBMAgamma0` <-
function(ensembleData, dates = NULL, trainingRule = list(length=NA,lag=NA),
         control = controlBMAgamma0(), warmStart = FALSE, 
         exchangeable = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 call <- match.call()

 trainingRule <- as.list(trainingRule)

 if (is.null(trainingRule$length) || is.na(trainingRule$length) ||
     is.null(trainingRule$lag) || is.na(trainingRule$lag))
   stop("length and lag must be specified in training rule")

 if (trainingRule$length < 0 || trainingRule$lag < 0) 
   stop("improperly specified training rule")

 ensMemNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable, ensembleGroups(ensembleData),
                                  nForecasts)

# remove instances missing all forecasts, obs or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 nObs <- ensembleNobs(ensembleData)
 if (!nObs) stop("no observations")

 ensDates <- ensembleDates(ensembleData)
 if (is.null(ensDates)) stop("dates unavailable")

 Dates <- as.character(ensDates)
 DATES <- sort(unique(Dates))

 if (trainingRule$length+trainingRule$lag > length(DATES)) 
   stop("insufficient training data")

 julianDATES <- ymdhTOjul(DATES)
 origin <- attr( julianDATES, "origin")
 incr <- min(1,min(diff(julianDATES))) ## incr may be fractional for hours

## dates that can be modeled by the training data (ignoring gaps)

 Jdates <- seq(from = julianDATES[trainingRule$length]+trainingRule$lag*incr,
               to = max(julianDATES)+trainingRule$lag*incr, by = incr)

## determine the modeling dates

 dl <- nchar(DATES[1])

 if (nullDates <- is.null(dates)) {

   dates <- julTOymdh(Jdates, origin = origin, dropHour = (dl == 8))

 }
 else {

   dates <- sort(unique(as.character(dates)))

   if (!all(dateCheck(dates))) 
     stop("improperly specified date(s) in dates argument")

   if (nchar(dates[1]) != dl)
     stop("format of dates argument does not match date format in data")

   if (any(dates < julTOymdh(min(Jdates),origin=origin,dropHour=(dl == 8)))) {
     stop("some dates precede the first training period")
   }

   if (any(dates > julTOymdh(max(Jdates),origin=origin,dropHour=(dl == 8)))) {
     warning("there are dates beyond the last training period")
   }

 }

 juliandates <- ymdhTOjul( dates, origin = origin)

 nDates <- length(dates)

 prob0coefs <- array( NA, c(3, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, dates))
 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, dates))
 varCoefs <- array( NA, c(2, nDates), dimnames = dates)
 weights <- array( NA, c(nForecasts, nDates),
                      dimnames = list(ensMemNames, dates))

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- rep(0, nDates)
 names(nIter) <- dates

 L <- length(juliandates)
 twin <- 1:trainingRule$length

## temp <- data.frame(julian = julianDATES,date = DATES)
## print(temp)

 cat("\n")
 l <- 0
 for(i in seq(along = juliandates)) {

    I <- (juliandates[i]-trainingRule$lag*incr) >= julianDATES
    if (!any(I)) stop("insufficient training data")

    j <- which(I)[sum(I)]

    if (j != l) {

      twin <- (j+1) - (1:trainingRule$length)

      D <- as.logical(match(Dates, DATES[twin], nomatch=0))
      if (!any(D)) stop("this should not happen")

      cat("modeling for date", dates[i], "...")

      fit <- fitBMAgamma0(ensembleData[D,], control = control,
                          exchangeable = exchangeable)

      l <- j ## last model fit
      trainTable[i] <- length(unique(Dates[D]))
      nIter[i] <- fit$nIter 
      if (warmStart) control$start$weights <- weights[,i]
      cat("\n")
     }
   else {
     trainTable[i] <- -trainTable[i-1]
     nIter[i] <- -nIter[i-1]
   }

   prob0coefs[,,i] <- fit$prob0coefs
   biasCoefs[,,i] <- fit$biasCoefs
   varCoefs[,i] <- fit$varCoefs
   weights[,i] <- fit$weights

 }

 structure(list(trainingRule = trainingRule, 
                prob0coefs = prob0coefs, biasCoefs = biasCoefs, 
                varCoefs = varCoefs, weights = weights, nIter = nIter,
                exchangeable = exchangeable,
                transformation = fit$transformation, 
                inverseTransformation = fit$inverseTransformation),
                call = match.call(), class = "ensembleBMAgamma0")
}
