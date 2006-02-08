"crpsANDmae.ensembleBMAgamma0" <-
function(object, ensembleData, nSamples = 10000, seed = NULL, popData = NULL, ...) 
{

 sampleFromDist <- function( x, n=10000)
{
  p <- x[1]
  q <- x[2]
  shape <- x[3]
  rate <- x[4]

  c(rep(0, round(p*n)), rgamma(round(q*n), shape=shape, rate=rate))
} 

 inverseLogit <- function(x) {
# logit function safeguared against underflow and overflow
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

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 if (!is.null(seed)) set.seed(seed)

 Q <- as.vector(quantileForecastBMA( object, ensembleData))

 ensDates <- ensembleDates(ensembleData)
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
 crpsSim <- sampleMedian <- sampleMed <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsNames(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 k <- 0
 for (j in J) {
    k <- k + 1
 
    if (any(is.na(WEIGHTS <- object$weights[,k]))) next

    I <- which(as.logical(match(Dates, DATES[j], nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       VAR <- object$varCoefs[1,k] + object$varCoefs[2,k]*f

       fTrans <- sapply(f, object$transformation)

       MEAN <- apply(rbind(1, fTrans) * object$biasCoefs[,,k], 2, sum)

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       if (is.null(popData)) {
         POP <- sapply(apply(rbind( 1, fTrans, f==0) * object$popCoefs[,,k],
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         POP <- sapply(apply(rbind( 1, fTrans, popi) * object$popCoefs[,,k],
                                    2,sum), inverseLogit)
       }

       SAMPLES <- sample( 1:nForecasts, size = nSamples, 
                          replace = TRUE, prob = WEIGHTS)

       tab <- table(SAMPLES)

       tab <- table(unlist(apply( cbind(as.numeric(names(tab)), tab), 1,
              function(nj,POP) {
         sample(c(nj[1],0), size = nj[2], replace = TRUE,
                prob = c(1-POP[nj[1]],POP[nj[1]]))}, 
                POP = POP)))

       if (length(tab) > 1) {
          S <- apply( cbind(as.numeric(names(tab[-1])), tab[-1]), 1,
              function(nj,SHAPE,RATE) 
                  rgamma(nj[2], shape=SHAPE[nj[1]], rate=RATE[nj[1]]),
                                        SHAPE = SHAPE, RATE = RATE)
           
# model is fit to the cube root of the forecast

         S <- sapply(as.vector(unlist(S)),
                           object$inverseTransformation)

         SAMPLES <- c(rep(0, tab[1]), S)
       }
       else SAMPLES <- rep(0,tab[1])

##     p <- POP*WEIGHTS
##     q <- (1-POP)*WEIGHTS

##     SAMPLES <- apply(rbind(p=p,q=q,shape=SHAPE,rate=RATE), 2, 
##                         sampleFromDist, n = nSamples)

##       SAMPLES <- sapply(as.vector(unlist(SAMPLES)), 
##                         object$inverseTransformation)

#      sampleMean[i] <- mean(SAMPLES) 
       sampleMedian[i] <- median(SAMPLES) 
  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

 L <- which(!is.na(crpsSim))
 l <- length(L)

 crpsSim <- mean(crpsSim[L])

 crpsCli <- sapply(obs[L], function(x,Y) mean(abs(Y-x)), Y = obs[L])
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns <- mean(apply(abs(sweep(ensembleData[L,],MARGIN=2,FUN ="-",STATS=obs))
                   ,1,mean)
  - apply(ensembleData[L,], 2, function(z,Z) 
                apply(abs(sweep(Z, MARGIN = 2, FUN = "-", STATS = z)),1,sum),
                  Z = ensembleData[L,])/(2*(nForecasts*nForecasts)))

 maeCli <- mean(abs(obs[L] - median(obs[L])))
 maeEns <- mean(abs(obs[L] - apply(ensembleData[L,], 1, median)))
 maeBMA <- mean(abs(obs[L] - Q))
 maeSim <- mean(abs(obs[L] - sampleMedian[L]))

 MAT <- matrix( NA, 2, 4)
 dimnames(MAT) <- list(c("CRPS", "MAE"), 
                       c("climatology", "ensemble", "BMA", "simBMA"))

 MAT["CRPS", "climatology"] <- crpsCli
 MAT["CRPS", "ensemble"] <- crpsEns
 MAT["CRPS", "simBMA"] <- crpsSim

 MAT["MAE", "climatology"] <- maeCli
 MAT["MAE", "ensemble"] <- maeEns
 MAT["MAE", "BMA"] <- maeBMA
 MAT["MAE", "simBMA"] <- maeSim

 MAT
}

