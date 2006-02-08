"cdfBMA.ensembleBMAnormal" <-
function(fit, ensembleData, values, dates = NULL, ...) 
{

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 nObs <- ensembleNobs(ensembleData)

 if (!is.null(dates)) {
   K <- match(dates, names(fit$dateTable), nomatch=0)
   if (any(!K) || !length(K)) 
     stop("parameters not available for a specified date")
   dateTable <- fit$dateTable[K]
 }
 else {
   dateTable <- fit$dateTable
   K <- 1:length(dateTable)
  }

 if (is.null(ensDates <- ensembleDates(ensembleData))) {
   if (length(dateTable) > 1) stop("date ambiguity")
   Dates <- rep(1,nObs)
   dates <- DATES <- 1
   L <- 1:nObs
 }
 else {
   if (!is.null(dates)) {
     L <- as.logical(match(dates, as.character(ensDates), nomatch=0))
     if (all(!L) || !length(L)) 
       stop("specified dates incompatible with ensemble data")
   }
   Dates <- as.numeric(ensDates)
   DATES <- sort(unique(Dates))
   L <- as.logical(match(as.character(ensDates), names(dateTable), nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   dates <- sort(unique(Dates[L]))
 }

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData), as.character(values)) 

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0 
 for (j in J) {

    l <- l + 1
    k <- K[l]

    if (any(is.na(WEIGHTS <- fit$weights[,k]))) next
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd[,k] 
          }
          else rep(fit$sd[k], nForecasts)

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,k], 2, sum)

       CDF[i,] <- sapply( values, normalBMAcdf,
                          WEIGHTS = WEIGHTS, MEAN = MEAN, SD = SD) 

    }

 }

 CDF[ L, , drop = FALSE]
}

