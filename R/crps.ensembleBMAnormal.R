`crps.ensembleBMAnormal` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 weps <- 1.e-4

 erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

 absExp <- function(mu, sig) 
  {
   (sqrt(2)* sig)*exp(-(mu/sig)^2/2)/sqrt(pi) + 
                       mu * erf((sqrt(2)*mu)/(2*sig))
  }

 matchITandFH(fit,ensembleData)

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

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    WEIGHTS <- fit$weights[,k]
    if (all(Wmiss <- is.na(WEIGHTS))) next
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd[,k] 
          } 
          else {
            rep(fit$sd[k], nForecasts)
          }

    VAR <- SD*SD

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,k], 2, sum)

       M <- is.na(f) | Wmiss

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W[!M] <- W[!M] / sum(W[!M])
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

         MEAN <- MEAN[!M]
         SD <- SD[!M]

         if (sum(!M) > 1) {
           SAMPLES <- sample((1:nForecasts)[!M],size=nSamples,
                             replace=TRUE,prob=W[!M]) 
         }
         else {
           SAMPLES <- rep( (1:nForecasts)[!M], nSamples)
         }

         tab <- table(SAMPLES)

         Z <- tab == 0

         tab <- tab[!Z]

         MEAN <- MEAN[!Z]
         SD <- SD[!Z]

         SAMPLES <- apply(cbind(seq(along = tab), tab), 1,
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

 crpsBMA <-  if (is.null(nSamples)) CRPS else crpsSim

 crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
 crpsCli <- crpsCli - mean(crpsCli)/2

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

 crpsEns <- crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts))

#cbind(climatology = crpsCli, ensemble = crpsEns, BMA = crpsBMA)
 cbind(ensemble = crpsEns, BMA = crpsBMA)
}

