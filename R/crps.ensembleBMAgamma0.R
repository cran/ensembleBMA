`crps.ensembleBMAgamma0` <-
function(fit, ensembleData, nSamples=10000, seed=NULL, dates=NULL, ...) 
{

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 if (is.null(nSamples)) nSamples <- 10000

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]
## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]
 
## match specified dates with dateTable in fit

 dateTable <- dimnames(fit$weights)[[2]]

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable)) 
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K)) 
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L)) 
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 obs <- ensembleVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- rep(NA, nObs)
 names(crpsSim) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    WEIGHTS <- fit$weights[,k]
    if (all(Wmiss <- is.na(WEIGHTS))) next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       VAR <- fit$varCoefs[1,k] + fit$varCoefs[2,k]*f

       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,,k], 2, sum)

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       PROB0 <- sapply(apply(rbind( 1, fTrans, f==0)*fit$prob0coefs[,,k],
                              2,sum), inverseLogit)
       RATE  <- RATE[!M]
       SHAPE <- SHAPE[!M]
       PROB0 <- PROB0[!M]

       if (sum(!M) > 1) {
         W <- WEIGHTS + weps
         W <- W[!M] / sum(W[!M])
         SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples, 
                           replace = TRUE, prob = W)
       }
       else {
         SAMPLES <- rep( (1:nForecasts)[!M], nSamples)
       }

       tab <- numeric(sum(!M))
       names(tab) <- (1:nForecasts)[!M]
       tabSamples <- table(SAMPLES)
       tab[names(tabSamples)] <- tabSamples

       Z <- tab == 0

       tab <- tab[!Z]

       RATE  <- RATE[!Z]
       SHAPE <- SHAPE[!Z]
       PROB0 <- PROB0[!Z]

       tab0 <- table(unlist(apply( cbind( seq(along = tab), tab), 1,
              function(nj,PROB0) {
         sample(c(nj[1],0), size = nj[2], replace = TRUE,
                prob = c(1-PROB0[nj[1]],PROB0[nj[1]]))}, 
                PROB0 = PROB0)))

       I <- as.numeric(names(tab0[-1]))
       tab[] <- 0
       tab[I] <- tab0[-1]

       Z <- tab == 0

       tab <- tab[!Z]

       if (length(tab)) {

          RATE  <- RATE[!Z]
          SHAPE <- SHAPE[!Z]

          S <- apply( cbind( seq(along = tab), tab), 1,
              function(nj,SHAPE,RATE) 
                  rgamma(nj[2], shape=SHAPE[nj[1]], rate=RATE[nj[1]]),
                                      SHAPE = SHAPE, RATE = RATE)
           
# model is fit to the cube root of the forecast

         S <- sapply(as.vector(unlist(S)), powinv, power = fit$power)

         SAMPLES <- c(rep(0, tab0[1]), S)
       }
       else SAMPLES <- rep(0,tab0[1])
  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

 crpsSim <- mean(crpsSim, na.rm = TRUE)

 crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
 crpsCli <- mean(crpsCli - mean(crpsCli)/2)

 crpsEns1 <- apply(abs(sweep(ensembleData,MARGIN=1,FUN ="-",STATS=obs))
                   ,1,mean,na.rm=TRUE)

 if (nrow(ensembleData) > 1) {
   crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
     apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
                  Z = ensembleData),1,sum, na.rm = TRUE)
 }
 else {
   crpsEns2 <- sum(sapply(as.vector(ensembleData), 
                   function(z,Z) sum( Z-z, na.rm = TRUE),
                   Z = as.vector(ensembleData)), na.rm = TRUE)
 }

 crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))

##c(climatology = crpsCli, ensemble = crpsEns, BMA = crpsSim)
c(ensemble = crpsEns, BMA = crpsSim)
}

