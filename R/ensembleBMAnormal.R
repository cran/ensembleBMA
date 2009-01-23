`ensembleBMAnormal` <-
function(ensembleData, trainingDays, dates = NULL, 
         control = controlBMAnormal(), exchangeable = NULL, minCRPS = FALSE)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 call <- match.call()

 warmStart <- FALSE

 if (is.list(trainingDays)) trainingsDays <- trainingsDays[[1]]

 if (length(trainingDays) > 1 || trainingDays <= 0 
       || (trainingDays - trunc(trainingDays)) != 0) 
   stop("trainingDays improperly specified")

 forecastHour <- ensembleFhour(ensembleData)
 lag <- ceiling( forecastHour / 24 )

 ensMemNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable, ensembleGroups(ensembleData),
                                  nForecasts)

# remove instances missing all forecasts, obs or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleValidDates(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 nObs <- ensembleNobs(ensembleData)
 if (!nObs) stop("no observations")

 ensDates <- ensembleValidDates(ensembleData)
 if (is.null(ensDates)) stop("dates unavailable")

 Dates <- as.character(ensDates)
 DATES <- sort(unique(Dates))

 if (trainingDays+lag > length(DATES)) 
   stop("insufficient training data")

 julianDATES <- ymdhTOjul(DATES)
 origin <- attr( julianDATES, "origin")
 incr <- min(1,min(diff(julianDATES))) ## incr may be fractional for hours

## dates that can be modeled by the training data (ignoring gaps)

 Jdates <- seq(from = julianDATES[trainingDays]+lag*incr,
               to = max(julianDATES)+lag*incr, by = incr)

## determine the modeling dates

 DATEShh <- getHH(DATES)

 if (length(DATEShh) != 1) stop("forecast hour in data should be unique")

 lD <- nchar(DATES[1])

 if (!(lD <- unique(sapply(DATES,nchar)))) 
    stop("all dates in data should have same character length")

 if (nullDates <- is.null(dates)) {

   dates <- julTOymdh(Jdates, origin = origin, dropHour = (lD == 8))

 }
 else {

   dates <- sort(unique(as.character(dates)))

   if (!all(dateCheck(dates))) 
     stop("improperly specified date(s) in dates argument")

   datesHH <- getHH(dates)

   if (length(datesHH) != 1) stop("forecast hour in dates should be unique")
   
   if (datesHH != DATEShh) stop("specified dates incompatible with data")

   if (!(ld <- unique(sapply(dates,nchar)))) 
     stop("all specified dates should have same character length")

   if (ld < lD) {
     dates <- sapply( dates, function(s) paste(s, "00", sep =""))
   }
   else if (ld < lD) {
     dates <- sapply( dates, function(s) substring(s, 1, 8))
   }

   if (any(dates < julTOymdh(min(Jdates),origin=origin,dropHour=(lD == 8)))) {
     stop("some dates precede the first training period")
   }

   if (any(dates > julTOymdh(max(Jdates),origin=origin,dropHour=(lD == 8)))) {
     warning("there are dates beyond the last training period")
   }

 }

 juliandates <- ymdhTOjul( dates, origin = origin)

 nDates <- length(dates)

 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                     dimnames = list(NULL, ensMemNames, dates))

 if (control$equalVariance) {
   sd <- rep(NA, nDates) 
   names(sd) <- dates
 }
 else {
   sd <- array(NA,c(nForecasts,nDates), 
               dimnames=list(ensMemNames, dates))
 }

 weights <- array( NA, c(nForecasts, nDates))
 dimnames(weights) <- list(ensMemNames, dates)

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- loglikelihood <- rep(0, nDates)
 names(nIter) <- names(loglikelihood) <- dates

 L <- length(juliandates)
 twin <- 1:trainingDays

## temp <- data.frame(julian = julianDATES,date = DATES)
## print(temp)

 K <- 1:nForecasts

 cat("\n")

 l <- 0
 for(i in seq(along = juliandates)) {

    I <- (juliandates[i]-lag*incr) >= julianDATES
    if (!any(I)) stop("insufficient training data")

    j <- which(I)[sum(I)]

    if (j != l) {

      twin <- (j+1) - (1:trainingDays)
      D <- as.logical(match(Dates, DATES[twin], nomatch=0))
      if (!any(D)) stop("this should not happen")
      cat("modeling for date", dates[i], "...")

      kNA <- apply(ensembleForecasts(ensembleData[D,]), 2, 
                   function(x) all(is.na(x)))

      if (any(kNA)) {
        if (!is.null(x <- exchangeable)) x <- exchangeable[-K[kNA]]
        fit <- fitBMAnormal(ensembleData[D,-K[kNA]], control = control,
                            exchangeable = x)
      }
      else {
        fit <- fitBMAnormal(ensembleData[D,], control = control,
                            exchangeable = exchangeable)
      }
  
      if (minCRPS) {

        CRPSobjective <- function(SD) {
            if (any(kNA)) {
              crpsNormal(SD, fit$weights, fit$biasCoefs, ensembleData[D,-K[kNA]])
            }
            else {
              crpsNormal(SD, fit$weights, fit$biasCoefs, ensembleData[D,])
            }
         }

        if (control$equalVariance) {
          fit$sd <- optimize(CRPSobjective, interval = c(0, 6*fit$sd))$minimum
        }
        else {
          opt <- optim(fit$sd, CRPSobjective, method = "BFGS")
          if (!opt$convergence) {
            fit$sd <- opt$par
          }
          else {
            warning("CRPS could not be minimized")
          }
        }
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

   biasCoefs[,K[!kNA],i] <- fit$biasCoefs
   weights[K[!kNA],i] <- fit$weights

   if (control$equalVariance) {
     sd[i] <- fit$sd 
   }
   else {
     sd[K[!kNA],i] <- as.vector(fit$sd)
   }

 }

 structure(list(training = list(days=trainingDays,lag= lag,table=trainTable),
                biasCoefs = biasCoefs, sd = sd, weights = weights,
                nIter = nIter, exchangeable = exchangeable),
                forecastHour = forecastHour, 
                initializationTime = ensembleItime(ensembleData),
                call = match.call(), class = "ensembleBMAnormal")
}

