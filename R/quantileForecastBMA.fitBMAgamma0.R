"quantileForecastBMA.fitBMAgamma0" <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, popData = NULL, ...) 
{

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 inverseLogit <- function(x) {
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  1/(1+exp(-x))
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  if (x >= log(.Machine$double.eps)) {
                    x <- exp(x)
                    x/(1+x)
                  }
                  else exp(x)
                }
                else 0
             }
            }

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (dim(popData) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")
 
 L <- 1:nObs

 if (!is.null(dates) && !is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

if (!any(is.na(WEIGHTS <- fit$weights))) {

   for (i in L) {
    
       f <- ensembleData[i,]

       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f

       fTrans <- sapply(f, fit$transformation)

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f==0) * fit$prob0coefs,
                                2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi) * fit$prob0coefs,
                                2,sum), inverseLogit)
       }

       MEAN <- apply(rbind(1, fTrans)*fit$biasCoefs, 2, sum)

       Q[i,] <- sapply( quantiles, gamma0BMAquant, WEIGHTS=WEIGHTS,
                       PROB0=PROB0, MEAN=MEAN, VAR=VAR)
   }
 }

  apply(Q, 2, fit$inverseTransformation)
}

