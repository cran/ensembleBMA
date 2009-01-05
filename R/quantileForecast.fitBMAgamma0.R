`quantileForecast.fitBMAgamma0` <-
function(fit, ensembleData, quantiles = 0.5, dates=NULL, ...) 
{ 

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)
 
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)

# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights
 if (!all(Wmiss <- is.na(WEIGHTS))) {

   for (i in 1:nObs) {
    
       f <- ensembleData[i,]       

       M <- is.na(f) | Wmiss

       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f

       fTrans <- sapply( f, power, fit$power)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f==0) * fit$prob0coefs,
                             2,sum), inverseLogit)

       MEAN <- apply(rbind(1, fTrans)*fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply( quantiles, quantBMAgamma0, WEIGHTS=W,
                        MEAN=MEAN[!M], VAR=VAR[!M],PROB0=PROB0[!M])
   }
 }

  apply(Q, 2, powinv, power = fit$power)
}

