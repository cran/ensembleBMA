plot.ensembleBMAgamma0 <-
function(x, ensembleData, dates = NULL, ask = TRUE, ...) 
{

 par(ask = ask)

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 matchITandFH(x,ensembleData)

 exchangeable <- x$exchangeable

 M <- matchEnsembleMembers(x,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]
 
## match specified dates with dateTable in fit

 dateTable <- dimnames(x$weights)[[2]]

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable)) 
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K)) 
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 nForecasts <- ensembleSize(ensembleData)

 obs <- ensembleVerifObs(ensembleData)
 
 ensembleData <- ensembleForecasts(ensembleData)
 
 l <- 0
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    WEIGHTS <- x$weights[,k]
     
    if (all(Wmiss <- is.na(WEIGHTS))) next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       VAR <- x$varCoefs[1,k] + x$varCoefs[2,k]*f
        
       fTrans <- sapply(f, powfun, power = x$power)

       MEAN <- apply(rbind(1, fTrans) * x$biasCoefs[,,k], 2, sum)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*x$prob0coefs[,,k],
                              2,sum), inverseLogit)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

        plotBMAgamma0(WEIGHTS = W, MEAN = MEAN[!M], VAR = VAR[!M], 
                      PROB0 = PROB0[!M], obs = obs[i],
                      exchangeable = exchangeable)

    }

 }

 invisible()
}

