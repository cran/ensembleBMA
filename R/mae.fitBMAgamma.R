`mae.fitBMAgamma` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)
 
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

 nForecasts <- ensembleSize(ensembleData) 

 Q <- as.vector(quantileForecast( fit, ensembleData))

 sampleMedian <- sampleMean <- predictiveMean <- rep(NA, nObs)
 names(sampleMedian) <- names(sampleMean) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS))) {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2

       fTrans <- sapply( f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M] / sum(W[!M])
       }

       predictiveMean[i] <- sum(W * MEAN[!M])

    if (!is.null(nSamples)) {
 
       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples, 
                          replace = TRUE, prob = W)

       tab <- table(SAMPLES)

          S <- apply( cbind(as.numeric(names(tab)), tab), 1,
              function(nj,SHAPE,RATE) 
                  rgamma(nj[2], shape=SHAPE[nj[1]], rate=RATE[nj[1]]),
                                        SHAPE = SHAPE[!M], RATE = RATE[!M])

         SAMPLES <- sapply(as.vector(unlist(S)), powinv, power = fit$power)

       sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
    }
 }

}

## maeCli <- mean(abs(obs - median(obs)))
## maeEns <- mean(abs(obs - apply(ensembleData, 1, median)))

 maeCli <- mean(abs(obs - mean(obs)))
 maeEns <- mean(abs(obs - apply(ensembleData, 1, mean)))

 if (is.null(nSamples)) {
## maeBMA <- mean(abs(obs - predictiveMean))
   maeBMA <- mean(abs(obs - Q))
 }
 else {
## maeBMA <- maeSim <- mean(abs(obs - sampleMean))
   maeBMA <- maeSim <- mean(abs(obs - sampleMedian))
 }

c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)
c(ensemble = maeEns, BMA = maeBMA)
}

