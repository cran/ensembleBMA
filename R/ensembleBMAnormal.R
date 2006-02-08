"ensembleBMAnormal" <-
function(ensembleData, dates = NULL, trainingRule = list(length=30, lag=2), 
         control = controlBMAnormal(), warmStart = FALSE, minCRPS = FALSE)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 call <- match.call()

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

 if (any(!J)) stop("in sufficient training data")

 I <- as.logical(match(Dates, dates, nomatch = 0))

 dateTable <- table(Dates[I])
 if (is.factor(ensDates)) 
   names(dateTable) <- levels(ensDates)[as.numeric(names(dateTable))]

 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                     dimnames = list(NULL, ensMemNames, names(dateTable)))

 if (control$equalVariance) {
   sd <- rep(NA, nDates) 
   names(sd) <- names(dateTable)
 }
 else {
   sd <- array(NA,c(nForecasts,nDates), 
               dimnames=list(ensMemNames, names(dateTable)))
 }

 weights <- array( NA, c(nForecasts, nDates))
 dimnames(weights) <- list(ensMemNames, names(dateTable))

 trainingRule <- as.list(trainingRule)

 offset <- 1 - trainingRule$lag - (1:trainingRule$length)

 k <- 0
 for (j in J) {
    k <- k + 1
    if (any(j + offset < 1)) next
    D <- as.logical(match(Dates, DATES[j+offset], nomatch=0))
    if (!any(D)) 
      stop("specified dates do not match dates in ensembleData")
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
    biasCoefs[,,k] <- fit$biasCoefs
    weights[,k] <- fit$weights
    if (warmStart) control$start$weights <- fit$weights
    if (control$equalVariance) sd[k] <- fit$sd else sd[,k] <- fit$sd 
 }

 structure(list(dateTable = dateTable, trainingRule = trainingRule, 
                biasCoefs = biasCoefs, sd = sd, weights = weights),
                call = match.call(), class = "ensembleBMAnormal")
}

