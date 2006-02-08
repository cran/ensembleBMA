"ensembleBMAgamma0" <-
function(ensembleData, dates = NULL, trainingRule = list(length=30, lag=2), 
         control = controlBMAgamma0(), warmStart = FALSE, 
         exchangeable = NULL, popData = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 call <- match.call()

 trainingRule <- as.list(trainingRule)

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)

 ensMemNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable, ensembleGroups(ensembleData),
                                  nForecasts)

 ensDates <- ensembleDates(ensembleData)
 if (is.null(ensDates)) stop("dates unavailable")

 Dates <- as.character(ensDates)
 DATES <- sort(unique(Dates))

 if (is.null(dates)) {

   if (trainingRule$length+trainingRule$lag > length(DATES)) 
     stop("insufficient training data")

   dates <- DATES[(trainingRule$length+trainingRule$lag):length(DATES)]
 }
 else {
   datesInput <- dates
   dates <- unique(dates)
   if (is.character(dates)) {
     I <- sapply(dates, function(d,D) which(d == as.character(D))[1],
                 D = ensDates)
     dates <- sort(Dates[I])
   }
 }

 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)
 M <- match(DATES, dates, nomatch = 0)

 if (any(!J)) {
   print(dates[!J])
   stop("specified dates not matched in data")
 }

 if (all(!M)) {
   print(datesInput)
   stop("specified dates not matched in data")
 }

 I <- as.logical(match(Dates, dates, nomatch = 0))

 dateTable <- table(Dates[I])
# if (is.factor(ensDates)) 
#   names(dateTable) <- levels(ensDates)[as.numeric(names(dateTable))]

 prob0coefs <- array( NA, c(3, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, names(dateTable)))
 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, names(dateTable)))
 varCoefs <- array( NA, c(2, nDates), dimnames = list(NULL, names(dateTable)))
 weights <- array( NA, c(nForecasts, nDates),
                      dimnames = list(ensMemNames, names(dateTable)))

 nIter <- dateTable
 nIter[] <- 0

 trainingRule <- as.list(trainingRule)

 offset <- 1 - trainingRule$lag - (1:trainingRule$length)

 if (any(M <- apply(outer(J, offset, "+"), 1, function(x) any(x < 1)))){
   print(dates[which(M)])
   stop("insufficient training data")
 }

 cat("\n")

 k <- 0
 for (j in J) {
    k <- k + 1
    if (any(j + offset < 1)) next
    D <- as.logical(match(Dates, DATES[j + offset], nomatch=0))
    if (!any(D)) 
      stop("specified dates do not match dates in ensembleData")
    cat("modeling for date", DATES[j], "...")
    if (is.null(popData)) {
      fit <- fitBMAgamma0(ensembleData[D,], control = control,
                          exchangeable = exchangeable)
    }
    else {
      fit <- fitBMAgamma0(ensembleData[D,], control = control,
                          exchangeable = exchangeable,
                          popData = lapply(popData, function(x,D) x[D,], D=D))
    }
    prob0coefs[,,k] <- fit$prob0coefs
    biasCoefs[,,k] <- fit$biasCoefs
    varCoefs[,k] <- fit$varCoefs
    weights[,k] <- fit$weights
    nIter[k] <- fit$nIter
    if (warmStart) control$start$weights <- fit$weights
    cat("\n")
 }

 structure(list(dateTable = dateTable, trainingRule = trainingRule, 
                prob0coefs = prob0coefs, biasCoefs = biasCoefs, 
                varCoefs = varCoefs, weights = weights, nIter = nIter,
                exchangeable = attr(ensembleData,"exchangeable"),
                transformation = fit$transformation, 
                inverseTransformation = fit$inverseTransformation),
                call = match.call(), class = "ensembleBMAgamma0")
}

