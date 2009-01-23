`crps.fitBMAnormal` <-
function(fit, ensembleData, nSamples = NULL, seed=NULL, dates=NULL, ...) 
{
 weps <- 1.e-4

 if(!is.null(dates)) warning("dates ignored")

 erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

 absExp <- function(mu, sig) 
  {
   (sqrt(2)* sig)*exp(-(mu/sig)^2/2)/sqrt(pi) + 
                       mu * erf((sqrt(2)*mu)/(2*sig))
  }

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

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
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
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)

       M <- is.na(f) | Wmiss

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) 
  #   - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))
  # *sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.

     if (is.null(nSamples)) {      
       crps1 <- crps2 <- 0

  # Begin computing the first term in the CRPS formula.  
  # This is a double sum since it is over w(i)*w(j) for all i and j.

       for (f1 in (1:nForecasts)[!M]) 
         {
          for (f2 in (1:nForecasts)[!M]) 
             {
              tvar <- VAR[f1] + VAR[f2]  # total variance
             tsd <- sqrt(tvar)          # total standard deviation
              tmean <- MEAN[f1] - MEAN[f2]
              temp <- absExp(tmean,tsd)
              term <- (W[f1]*W[f2])*temp
              crps1 <- crps1 + term
             }
           tvar <- VAR[f1]              # total variance
           tsd <- sqrt(tvar)            # total standard deviation
           tmean <- MEAN[f1] - obs[i]
           crps2 <- crps2 + W[f1]*absExp(tmean,tsd)
        }

    # Using Szekely's expression for the CRPS, 
    # the first sum and second are put together to compute the CRPS.

       CRPS[i]  <- crps2 - crps1/2     
     }
     else {

       if (sum(!M) > 1) {
         SAMPLES <- sample((1:nForecasts)[!M],size=nSamples,
                           replace=TRUE,prob=W[!M]) 
       }
       else {
         SAMPLES <- rep( (1:nForecasts)[!M], nSamples)
       }

       tab <- table(SAMPLES)
       SAMPLES <- apply(cbind(as.numeric(names(tab)), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN[!M], SD = SD[!M])

       SAMPLES <- as.vector(unlist(SAMPLES))

#      sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1 <-  mean(abs(diff(sample(SAMPLES))))
       crps2  <- mean(abs(SAMPLES - obs[i])) 
       crpsSim[i]  <- crps2 - crps1/2
    }
  }
 }

 crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
 crpsCli <- crpsCli - mean(crpsCli)/2

 crpsEns1 <- apply(abs(sweep(ensembleData,MARGIN=1,FUN ="-",STATS=obs)),
                   1, mean, na.rm = TRUE)
 crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
       apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
                  Z = ensembleData), 1, sum, na.rm = TRUE)

 crpsEns <- crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts))

 crpsBMA <- if (is.null(nSamples)) CRPS else crpsSim

#cbind(climatology = crpsCli, ensemble = crpsEns, BMA = crpsBMA)
 cbind(ensemble = crpsEns, BMA = crpsBMA)
}

