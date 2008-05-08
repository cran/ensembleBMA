`cdf.fitBMAgamma0` <-
function(fit, ensembleData, values, dates=NULL, ...) 
{
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData), as.character(values)) 

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f
        
       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs,
                              2,sum), inverseLogit)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       CDF[i,] <- sapply( sapply( values, fit$transformation), 
                          cdfBMAgamma0, WEIGHTS = W, PROB0 = PROB0[!M], 
                          MEAN = MEAN[!M], VAR = VAR[!M]) 
    }

 }

 CDF
}

