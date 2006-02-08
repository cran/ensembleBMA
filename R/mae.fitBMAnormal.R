"mae.fitBMAnormal" <-
function(fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, ...) 
{
## contains CRPS code for historic reasons

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

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")

 L <- 1:nObs

 if (!is.null(dates) && !is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

 Q <- as.vector(quantileForecastBMA( fit, ensembleData))

 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- sampleMedian <- sampleMean <- predictiveMean <- rep(NA,nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)


 if (!any(is.na(WEIGHTS <- fit$weights))) {
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd
          } 
          else {
            rep(fit$sd, nForecasts)
          }

    VAR <- SD*SD

    for (i in L) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)

       predictiveMean[i] <- sum(WEIGHTS*MEAN)

 if (!is.null(nSamples)) {      

  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) 
  #   - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))
  # *sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.

##       crps1 <- crps2 <- 0

  # Begin computing the first term in the CRPS formula.  
  # This is a double sum since it is over w(i)*w(j) for all i and j.

##       for (f1 in 1:nForecasts) 
##         {
##          for (f2 in 1:nForecasts) 
##             {
##              tvar <- VAR[f1] + VAR[f2]  # total variance
##             tsd <- sqrt(tvar)          # total standard deviation
##              tmean <- MEAN[f1] - MEAN[f2]
##              temp <- absExp(tmean,tsd)
##              term <- (WEIGHTS[f1]*WEIGHTS[f2])*temp
##              crps1 <- crps1 + term
##             }
##           tvar <- VAR[f1]              # total variance
##           tsd <- sqrt(tvar)            # total standard deviation
##           tmean <- MEAN[f1] - obs[i]
##           crps2 <- crps2 + WEIGHTS[f1]*absExp(tmean,tsd)
##        }

    # Using Szekely's expression for the CRPS, 
    # the first sum and second are put together to compute the CRPS.

##       CRPS[i]  <- crps2 - crps1/2     


       SAMPLES <- sample(1:nForecasts,size=nSamples,replace=TRUE,prob=WEIGHTS) 
       tab <- table(SAMPLES)
       SAMPLES <- apply(cbind(as.numeric(names(tab)), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN, SD = SD)

       SAMPLES <- as.vector(unlist(SAMPLES))

       sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1 <-  mean(abs(diff(sample(SAMPLES))))
       crps2  <- mean(abs(SAMPLES - obs[i])) 
       crpsSim[i]  <- crps2 - crps1/2
    }
   }
 }

## L <- which(!is.na(crpsSim))
## l <- length(L)

## crpsSim <- mean(crpsSim[L])

 crpsCli <- sapply(obs[L], function(x,Y) mean(abs(Y-x)), Y = obs[L])
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns1 <- apply(abs(sweep(ensembleData[L,],MARGIN=1,FUN ="-",STATS=obs[L]))
                   ,1,mean)
 crpsEns2 <- apply(apply(ensembleData[L,], 2, function(z,Z) 
                apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum),
                  Z = ensembleData[L,]),1,sum)

 crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))

## maeCli <- mean(abs(obs[L] - median(obs[L])))
## maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, median)))

 maeCli <- mean(abs(obs[L] - mean(obs[L])))

 maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, mean)))

 if (is.null(nSamples)) {
## maeBMA <- mean(abs(obs[L] - Q[!is.na(Q)]))
## maeBMA <- mean(abs(obs[L] - Q))
   maeBMA <- mean(abs(obs[L] - predictiveMean[L]))
 } 
 else {
##  maeBMA <- maeSim <- mean(abs(obs[L] - sampleMedian[L]))
    maeBMA <- maeSim <- mean(abs(obs[L] - sampleMean[L]))
 }

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
##c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)
c(ensemble = maeEns, BMA = maeBMA)
}

