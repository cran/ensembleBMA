"crpsANDmae.ensembleBMAnormal" <-
function(object, ensembleData, nSamples = 10000, seed = NULL, ...) 
{

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

absExp <- function(mu, sig) 
  {
   (sqrt(2)* sig)*exp(-(mu/sig)^2/2)/sqrt(pi) + 
                       mu * erf((sqrt(2)*mu)/(2*sig))
  }

 if (!is.null(seed)) set.seed(seed)

 Q <- as.vector(quantileForecastBMA( object, ensembleData))

 ensDates <- ensembleData$dates
 Dates <- as.numeric(ensDates)
 DATES <- sort(unique(Dates))

 K <- sapply(names(object$dateTable), function(d,D) 
                        which(d == as.character(D))[1],
                            D = ensDates)
 dates <- sort(Dates[K])
 nDates <- length(dates)

 J <- match(dates, DATES, nomatch = 0)

 if (any(!J)) stop("specified dates not matched in data")

 nForecasts <- ensembleSize(ensembleData) 

 obs <- ensembleVerifObs(ensembleData)
 nObs <- length(obs)

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsNames(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 k <-0
 for (j in J) {
     k <- k + 1

    if (any(is.na(WEIGHTS <- object$weights[,k]))) next
     
    SD <- if (!is.null(dim(object$sd))) {
            object$sd[,k] 
          } 
          else {
            rep(object$sd[k], nForecasts)
          }

    VAR <- SD*SD

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * object$biasCoefs[,,k], 2, sum)

  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) 
  #   - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))
  # *sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.
      
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

       SAMPLES <- sample(1:nForecasts,size=nSamples,replace=TRUE,prob=WEIGHTS) 
       tab <- table(SAMPLES)
       SAMPLES <- apply(cbind(as.numeric(names(tab)), tab), 1,
                     function(nj,MEAN,SD) rnorm(nj[2],MEAN[nj[1]],SD[nj[1]]),
                        MEAN = MEAN, SD = SD)

       SAMPLES <- as.vector(unlist(SAMPLES))

#      sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1 <-  mean(abs(diff(sample(SAMPLES))))
       crps2  <- mean(abs(SAMPLES - obs[i])) 
       crpsSim[i]  <- crps2 - crps1/2
    }
 }

 L <- which(!is.na(crpsSim))
 l <- length(L)

 crpsSim <- mean(crpsSim[L])

 crpsCli <- sapply(obs[L], function(x,Y) mean(abs(Y-x)), Y = obs[L])
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns1 <- apply(abs(sweep(ensembleData[L,],MARGIN=1,FUN ="-",STATS=obs[L]))
                   ,1,mean)
 crpsEns2 <- apply(apply(ensembleData[L,], 2, function(z,Z) 
                apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum),
                  Z = ensembleData[L,]),1,sum)

 crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))

 maeSim <- mean(abs(obs[L] - sampleMedian[L]))

 maeCli <- mean(abs(obs[L] - median(obs[L])))

 maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, median)))
 maeBMA <- mean(abs(obs[L] - Q[!is.na(Q)]))

 MAT <- matrix( NA, 2, 4)
 dimnames(MAT) <- list(c("CRPS", "MAE"), 
                       c("climatology", "ensemble", "BMA", "simBMA"))

 MAT["CRPS", "climatology"] <- crpsCli
 MAT["CRPS", "ensemble"] <- crpsEns
 MAT["CRPS", "simBMA"] <- crpsSim
 MAT["CRPS", "BMA"] <- mean(CRPS,na.rm=TRUE)

 MAT["MAE", "climatology"] <- maeCli
 MAT["MAE", "ensemble"] <- maeEns
 MAT["MAE", "simBMA"] <- maeSim
 MAT["MAE", "BMA"] <- maeBMA

 MAT
}

