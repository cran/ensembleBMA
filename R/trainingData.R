"trainingData" <-
function(ensembleData, date, trainingRule = list(length = 30, lag = 2))
{

 if (length(date) > 1) stop("one day only")

 nObs <- ensembleVerifObs(ensembleData)
 ensNames <- ensembleMemberNames(ensembleData)
 nForecasts <- length(ensNames)

 ensDates <- ensembleDates(ensembleData)

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
 offset <- 1 - trainingRule$lag - (1:trainingRule$length)

 D <- as.logical(match(Dates, DATES[j+offset], nomatch=0))

 ensembleData[D,]
}

