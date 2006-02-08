"mae.ensembleBMAgamma0" <-
function(fit, ensembleData, dates=NULL, nSamples=10000, seed=NULL, 
         popData=NULL, ...) 
{
## contains CRPS code for historic reasons

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 sampleFromDist <- function( x, n=10000)
{
  p <- x[1]
  q <- x[2]
  shape <- x[3]
  rate <- x[4]

  c(rep(0, round(p*n)), rgamma(round(q*n), shape=shape, rate=rate))
} 

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 inverseLogit <- function(x) {
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  1/(1+exp(-x))
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  if (x >= log(.Machine$double.eps)) {
                    x <- exp(x)
                    x/(1+x)
                  }
                  else exp(x)
                }
                else 0
             }
            }

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
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

 crpsSim <- sampleMedian <- sampleMean <- predictiveMean <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 if (!is.null(nSamples)) {

 l <- 0
 for (j in J) {

    l <- l + 1
    k <- K[l]

    if (any(is.na(WEIGHTS <- fit$weights[,k]))) next

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       VAR <- fit$varCoefs[1,k] + fit$varCoefs[2,k]*f

       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,,k], 2, sum)

       predictiveMean[i] <- sum(WEIGHTS*MEAN)

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f==0) * fit$prob0coefs[,,k],
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi) * fit$prob0coefs[,,k],
                                    2,sum), inverseLogit)
       }

       SAMPLES <- sample( 1:nForecasts, size = nSamples, 
                          replace = TRUE, prob = WEIGHTS)

       tab <- table(SAMPLES)

       tab <- table(unlist(apply( cbind(as.numeric(names(tab)), tab), 1,
              function(nj,PROB0) {
         sample(c(nj[1],0), size = nj[2], replace = TRUE,
                prob = c(1-PROB0[nj[1]],PROB0[nj[1]]))}, 
                PROB0 = PROB0)))

       if (length(tab) > 1) {
          S <- apply( cbind(as.numeric(names(tab[-1])), tab[-1]), 1,
              function(nj,SHAPE,RATE) 
                  rgamma(nj[2], shape=SHAPE[nj[1]], rate=RATE[nj[1]]),
                                        SHAPE = SHAPE, RATE = RATE)
           
# model is fit to the cube root of the forecast

         S <- sapply(as.vector(unlist(S)),
                           fit$inverseTransformation)

         SAMPLES <- c(rep(0, tab[1]), S)
       }
       else SAMPLES <- rep(0,tab[1])

##     p <- PROB0*WEIGHTS
##     q <- (1-PROB0)*WEIGHTS

##     SAMPLES <- apply(rbind(p=p,q=q,shape=SHAPE,rate=RATE), 2, 
##                         sampleFromDist, n = nSamples)

##       SAMPLES <- sapply(as.vector(unlist(SAMPLES)), 
##                         fit$inverseTransformation)

#      sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
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
   maeBMA <- mean(abs(obs[L] - Q))
 }
 else {
   maeBMA <-  maeSim <- mean(abs(obs[L] - sampleMedian[L]))
 }

## MAT <- matrix( NA, 2, 4)
## dimnames(MAT) <- list(c("CRPS", "MAE"), 
##                       c("climatology", "ensemble", "BMA", "simBMA"))

## MAT["CRPS", "climatology"] <- crpsCli
## MAT["CRPS", "ensemble"] <- crpsEns
## MAT["CRPS", "simBMA"] <- crpsSim

## MAT["MAE", "climatology"] <- maeCli
## MAT["MAE", "ensemble"] <- maeEns
## MAT["MAE", "BMA"] <- maeBMA
## MAT["MAE", "simBMA"] <- maeSim

## c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)
 c(ensemble = maeEns, BMA = maeBMA)
}

