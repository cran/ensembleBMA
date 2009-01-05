`pit.fitBMAgamma` <-
function(fit, ensembleData, dates=NULL, ...) 
{
 
 powfun <- function(x,power) x^power

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 PIT <- numeric(nObs)
 names(PIT) <- ensembleObsLabels(ensembleData)

 obs <- sapply( ensembleVerifObs(ensembleData), powfun, power = fit$power)
 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2
        
       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       PIT[i] <- cdfBMAgamma( obs[i], WEIGHTS = W, 
                          MEAN = MEAN[!M], VAR = VAR[!M]) 
    }

 }

 PIT
}

