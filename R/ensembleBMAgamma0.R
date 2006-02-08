"ensembleBMAgamma0" <-
function(ensembleData, dates = NULL, trainingRule = list(length=30, lag=2), 
         control = controlBMAgamma0(), warmStart = FALSE, popData = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 call <- match.call()

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)
 ensMemNames <- ensembleMemberNames(ensembleData)
 nForecasts <- length(ensMemNames)

 ensDates <- ensembleDates(ensembleData)
 Dates <- as.numeric(ensDates)
 DATES <- sort(unique(Dates))

 if (is.null(dates)) {
   dates <- DATES[(trainingRule$length+trainingRule$lag):length(DATES)]
 }
 else {
   dates <- unique(dates)
   if (is.character(dates)) {
     I <- sapply(dates, function(d,D) which(d == as.character(D))[1],
                 D = ensDates)
     dates <- sort(Dates[I])
   }
 }

 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("insufficient training data")

 I <- as.logical(match(Dates, dates, nomatch = 0))

 dateTable <- table(Dates[I])
 if (is.factor(ensDates)) 
   names(dateTable) <- levels(ensDates)[as.numeric(names(dateTable))]

 popCoefs <- array( NA, c(3, nForecasts, nDates),
                      dimnames = list(NULL, NULL, names(dateTable)))
 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                      dimnames = list(NULL, ensMemNames, names(dateTable)))
 varCoefs <- array( NA, c(2, nDates), dimnames = list(NULL, names(dateTable)))
 weights <- array( NA, c(nForecasts, nDates),
                      dimnames = list(ensMemNames, names(dateTable)))

 trainingRule <- as.list(trainingRule)
 offset <- 1 - trainingRule$lag - (1:trainingRule$length)

 k <- 0
 for (j in J) {
    k <- k + 1
    if (any(j + offset < 1)) next
    D <- as.logical(match(Dates, DATES[j + offset], nomatch=0))
    if (!any(D)) 
      stop("specified dates do not match dates in ensembleData")
    if (is.null(popData)) {
      fit <- fitBMAgamma0(ensembleData[D,], control = control)
    }
    else {
      fit <- fitBMAgamma0(ensembleData[D,], control = control,
                 popData = lapply(popData, function(x,D) x[D,], D=D))
    }
    popCoefs[,,k] <- fit$popCoefs
    biasCoefs[,,k] <- fit$biasCoefs
    varCoefs[,k] <- fit$varCoefs
    weights[,k] <- fit$weights
    if (warmStart) control$start$weights <- fit$weights
 }

 structure(list(dateTable = dateTable, trainingRule = trainingRule, 
      popCoefs = popCoefs, biasCoefs = biasCoefs, varCoefs = varCoefs,
      weights = weights, transformation = fit$transformation, 
      inverseTransformation = fit$inverseTransformation),
      call = match.call(), class = "ensembleBMAgamma0")
}

