`mae.fitBMAnormal` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts or obs

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

 nObs <- ensembleNobs(ensembleData)

 if (!is.null(seed)) set.seed(seed)

 Q <- as.vector(quantileForecast( fit, ensembleData))

 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- sampleMedian <- sampleMean <- predictiveMean <- rep(NA,nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights
     
 if (!all(Wmiss <- is.na(WEIGHTS))) {
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd
          } 
          else {
            rep(fit$sd, nForecasts)
          }

    VAR <- SD*SD

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M] / sum(W[!M])
       }
       predictiveMean[i] <- sum(W*MEAN[!M])

       if (!is.null(nSamples)) {      

         SAMPLES <- sample((1:nForecasts)[!M],size=nSamples,
                           replace=TRUE,prob=W) 
  
         tab <- table(SAMPLES)
         SAMPLES <- apply(cbind(as.numeric(names(tab)), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN[!M], SD = SD[!M])

         SAMPLES <- as.vector(unlist(SAMPLES))

         sampleMean[i] <- mean(SAMPLES) 
         sampleMedian[i] <- median(SAMPLES) 
  
       }
     }
 }

## maeCli <- mean(abs(obs - median(obs)))
## maeEns <- mean(abs(obs - apply(ensembleData, 1, median)))

 maeCli <- mean(abs(obs - mean(obs)))

 maeEns <- mean(abs(obs - apply(ensembleData, 1, mean, na.rm = TRUE)))

 if (is.null(nSamples)) {
##  maeBMA <- mean(abs(obs - predictiveMean))
    maeBMA <- mean(abs(obs - Q))
 } 
 else {
##  maeBMA <- maeSim <- mean(abs(obs - sampleMean))
    maeBMA <- maeSim <- mean(abs(obs - sampleMedian))
 }

c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)
c(ensemble = maeEns, BMA = maeBMA)
}

