`trainingData` <-
function(ensembleData, date = NULL, trainingRule = list(length=NA, lag=NA), 
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

 if (nullDates <- is.null(dates)) stop("date must be specified")

 dates <- date


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

 juliandates <- ymdhTOjul( dates, origin = origin)

 nDates <- length(dates)

 L <- length(juliandates)
 twin <- 1:trainingRule$length

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
  
     DATA <- ensembleData[D,]

    }
 }

 attr( DATA, "exchangeable") <- exchangeable
 DATA
}

