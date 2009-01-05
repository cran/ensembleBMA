`crps.fitBMAgamma` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 if (is.null(nSamples)) nSamples <- 10000

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts or obs

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]

 if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

#nObs <- length(obs) 
 nObs <- ensembleNobs(ensembleData)

 if (!is.null(seed)) set.seed(seed)

 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- rep(NA, nObs)
 names(crpsSim) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS))) {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2

       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       if (sum(!M) > 1 ) {  
         W <- WEIGHTS + weps
         W <- W[!M]/sum(W[!M])
         SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples, 
                            replace = TRUE, prob = W)
       }
       else {
         SAMPLES <- rep( (1:nForecasts)[!M], nSamples)
       }

       tab <- table(SAMPLES)

       S <- apply( cbind(as.numeric(names(tab)), tab), 1,
              function(nj,SHAPE,RATE) 
                  rgamma(nj[2], shape=SHAPE[nj[1]], rate=RATE[nj[1]]),
                                        SHAPE = SHAPE[!M], RATE = RATE[!M])

       SAMPLES <- sapply(as.vector(unlist(S)), powinv, power = fit$power)

  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

 crpsSim <- mean(crpsSim)

 crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns1 <- apply(abs(sweep(ensembleData,MARGIN=1,FUN ="-",STATS=obs)),
                   1, mean, na.rm = TRUE)
 crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
       apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
                  Z = ensembleData), 1, sum, na.rm = TRUE)
 crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))

##c(climatology = crpsCli, ensemble = crpsEns, BMA = crpsSim)
c(ensemble = crpsEns, BMA = crpsSim)
}

