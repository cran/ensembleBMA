"crps.ensembleBMAnormal" <-
function(fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, ...) 
{
##
## has MAE also for historical reasons
##

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

absExp <- function(mu, sig) 
  {
   (sqrt(2)* sig)*exp(-(mu/sig)^2/2)/sqrt(pi) + 
                       mu * erf((sqrt(2)*mu)/(2*sig))
  }

 if (!is.null(seed)) set.seed(seed)

 if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

 nObs <- length(obs)

 Q <- as.vector(quantileForecastBMA( fit, ensembleData, dates = dates))

 if (!is.null(dates)) {
   K <- match(dates, names(fit$dateTable), nomatch=0)
   if (any(!K) || !length(K)) 
     stop("parameters not available for a specified date")
   dateTable <- fit$dateTable[K]
 }
 else {
   dateTable <- fit$dateTable
   K <- 1:length(dateTable)
 }

 if (is.null(ensDates <- ensembleDates(ensembleData))) {
   if (length(dateTable) > 1) stop("date ambiguity")
   Dates <- rep(1,nObs)
   dates <- DATES <- 1
   L <- 1:nObs
 }
 else {
   if (!is.null(dates)) {
     L <- as.logical(match(dates, as.character(ensDates), nomatch=0))
     if (all(!L) || !length(L)) 
       stop("specified dates incompatible with ensemble data")
   }
   Dates <- as.numeric(ensDates)
   DATES <- sort(unique(Dates))
   L <- as.logical(match(as.character(ensDates), names(dateTable), nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   dates <- sort(unique(Dates[L]))
 }

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 nForecasts <- ensembleSize(ensembleData) 

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (j in J) {
    l <- l + 1
    k <- K[l]

    if (any(is.na(WEIGHTS <- fit$weights[,k]))) next
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd[,k] 
          } 
          else {
            rep(fit$sd[k], nForecasts)
          }

    VAR <- SD*SD

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,k], 2, sum)

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

         for (f1 in 1:nForecasts) 
            {
            for (f2 in 1:nForecasts) 
               {
                tvar <- VAR[f1] + VAR[f2]  # total variance
                tsd <- sqrt(tvar)          # total standard deviation
                tmean <- MEAN[f1] - MEAN[f2]
                temp <- absExp(tmean,tsd)
                term <- (WEIGHTS[f1]*WEIGHTS[f2])*temp
                crps1 <- crps1 + term
               }
            tvar <- VAR[f1]              # total variance
            tsd <- sqrt(tvar)            # total standard deviation
            tmean <- MEAN[f1] - obs[i]
            crps2 <- crps2 + WEIGHTS[f1]*absExp(tmean,tsd)
         }

    # Using Szekely's expression for the CRPS, 
    # the first sum and second are put together to compute the CRPS.

         CRPS[i]  <- crps2 - crps1/2     
        }
       else {
         SAMPLES <- sample(1:nForecasts,size=nSamples,
                           replace=TRUE,prob=WEIGHTS) 
         tab <- table(SAMPLES)
         SAMPLES <- apply(cbind(as.numeric(names(tab)), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN, SD = SD)

         SAMPLES <- as.vector(unlist(SAMPLES))

#        sampleMean[i] <- mean(SAMPLES) 
         sampleMedian[i] <- median(SAMPLES) 
  
         crps1 <-  mean(abs(diff(sample(SAMPLES))))
         crps2  <- mean(abs(SAMPLES - obs[i])) 
         crpsSim[i]  <- crps2 - crps1/2
       }
    }
 }

## L <- which(!is.na(crpsSim))
## l <- length(L)
 crpsCli <- sapply(obs[L], function(x,Y) mean(abs(Y-x)), Y = obs[L])
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns1 <- apply(abs(sweep(ensembleData[L,],MARGIN=1,FUN ="-",STATS=obs[L]))
                   ,1,mean)
 crpsEns2 <- apply(apply(ensembleData[L,], 2, function(z,Z) 
                apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum),
                  Z = ensembleData[L,]),1,sum)

 crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))

 maeCli <- mean(abs(obs[L] - median(obs[L])))

 maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, median)))
 maeBMA <- mean(abs(obs[L] - Q[!is.na(Q)]))

## MAT <- matrix( NA, 2, 4)
## dimnames(MAT) <- list(c("CRPS", "MAE"), 
##                       c("climatology", "ensemble", "BMA", "simBMA"))
## MAT["CRPS", "climatology"] <- crpsCli
## MAT["CRPS", "ensemble"] <- crpsEns
## MAT["CRPS", "simBMA"] <- crpsSim
## MAT["CRPS", "BMA"] <- mean(CRPS,na.rm=TRUE)

## MAT["MAE", "climatology"] <- maeCli
## MAT["MAE", "ensemble"] <- maeEns
## MAT["MAE", "simBMA"] <- maeSim
## MAT["MAE", "BMA"] <- maeBMA

## MAT

 crpsBMA <- if (is.null(nSamples)) mean(CRPS,na.rm=TRUE) else mean(crpsSim[L])

 c(climatology = crpsCli, ensemble = crpsEns, BMA = crpsBMA)
 c(ensemble = crpsEns, BMA = crpsBMA)
}

