`trainingData` <-
function(ensembleData, date, trainingRule = list(length = NA, lag = NA))
{

 if (length(date) > 1) stop("one day only")

# remove instances missing all forecasts, obs, or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]

#nObs <- length(y)
 nObs <- ensembleNobs(ensembleData)
 if (!nObs) stop("no observations")

 ensNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensNames)

 ensDates <- ensembleDates(ensembleData)
 if (is.null(ensDates)) stop("dates unavailable")

 if (is.numeric(date)) {
   Dates <- as.numeric(ensDates)
 }
 else {
   Dates <- as.character(ensDates)
 }

 DATES <- unique(Dates)

 j <- match(date, DATES, nomatch = 0)

 if (any(!j)) 
   stop("specified date does not match dates in ensembleData")

 trainingRule <- as.list(trainingRule)
 if (is.null(trainingRule$length) || is.na(trainingRule$length) ||
     is.null(trainingRule$lag) || is.na(trainingRule$lag))
   stop("length and lag must be specified in training rule")

 offset <- 1 - trainingRule$lag - (1:trainingRule$length)

 if (any(j + offset < 1)) stop("insufficient training data")
 D <- as.logical(match(Dates, DATES[j+offset], nomatch=0))

 ensembleData[D,]
}

