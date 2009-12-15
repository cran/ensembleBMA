plot.fitBMAnormal <-
function(x, ensembleData, dates=NULL, ...) 
{

 exchangeable <- x$exchangeable

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(x,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]
 
# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]
 
 nObs <- nrow(ensembleData)

 obs <- ensembleVerifObs(ensembleData)

 nForecasts <- ensembleSize(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- x$weights
     
 if (!all(Wmiss <- is.na(WEIGHTS))) {
     
    SD <- if (!is.null(dim(x$sd))) {
            x$sd
          }
          else rep(x$sd, nForecasts)

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * x$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       plotBMAnormal( WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M],
                      obs = obs[i], exchangeable = exchangeable)

    }

 }

 invisible()
}

