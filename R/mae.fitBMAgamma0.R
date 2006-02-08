"mae.fitBMAgamma0" <-
function(fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, 
         popData=NULL, ...) 
{
## contains CRPS code for historical reasons

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

 crpsSim <- sampleMedian <- sampleMean <- predictiveMean <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 if (!any(is.na(WEIGHTS <- fit$weights))) {

    for (i in L) {
    
       f <- ensembleData[i,]
     
       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f

       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       predictiveMean[i] <- sum(WEIGHTS * MEAN)

 if (!is.null(nSamples)) {
 
       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f==0) * fit$prob0coefs,
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi) * fit$prob0coefs,
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

##     p <- POP*WEIGHTS
##     q <- (1-POP)*WEIGHTS

##     SAMPLES <- apply(rbind(p=p,q=q,shape=SHAPE,rate=RATE), 2, 
##                         sampleFromDist, n = nSamples)

##       SAMPLES <- sapply(as.vector(unlist(SAMPLES)), 
##                         fit$inverseTransformation)

       sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

 L <- which(!is.na(crpsSim))
 l <- length(L)

 crpsSim <- mean(crpsSim[L])
}

 crpsCli <- sapply(obs[L], function(x,Y) mean(abs(Y-x)), Y = obs[L])
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns <- mean(apply(abs(sweep(ensembleData[L,],MARGIN=2,FUN ="-",STATS=obs))
                   ,1,mean)
  - apply(ensembleData[L,], 2, function(z,Z) 
                apply(abs(sweep(Z, MARGIN = 2, FUN = "-", STATS = z)),1,sum),
                  Z = ensembleData[L,])/(2*(nForecasts*nForecasts)))

## maeCli <- mean(abs(obs[L] - median(obs[L])))
## maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, median)))

 maeCli <- mean(abs(obs[L] - mean(obs[L])))
 maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, mean)))

 if (is.null(nSamples)) {
## maeBMA <- mean(abs(obs[L] - predictiveMean[L]))
   maeBMA <- mean(abs(obs[L] - Q))
 }
 else {
## maeBMA <- maeSim <- mean(abs(obs[L] - sampleMean[L]))
   maeBMA <- maeSim <- mean(abs(obs[L] - sampleMedian[L]))
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

## MAT

##c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)
c(ensemble = maeEns, BMA = maeBMA)
}

