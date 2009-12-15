crps.ensembleBMAnormal <-
function(fit, ensembleData, nSamples=10000, seed=NULL, dates=NULL, ...) 
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

 Q <- as.vector(quantileForecast( fit, ensembleData, dates = dates))
 if (any(is.na(Q))) stop("NAs in forecast") # fix like ensembleBMAgamma0

 obs <- ensembleVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData) 
 members <- ensembleMemberLabels(ensembleData)

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

  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) 
  #   - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))
  # *sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.

       if (is.null(nSamples)) {      

         W <- WEIGHTS
         if (any(M)) {
           W <- W + weps
           W[!M] <- W[!M] / sum(W[!M])
         }

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
                 crps2 <- crps2 + term
                }
              tvar <- VAR[f1]              # total variance
              tsd <- sqrt(tvar)            # total standard deviation
              tmean <- MEAN[f1] - obs[i]
              crps1 <- crps1 + W[f1]*absExp(tmean,tsd)
            }

    # Using Szekely's expression for the CRPS, 
    # the first sum and second are put together to compute the CRPS.

         CRPS[i]  <- crps1 - crps2/2     
        }
      else {

          W <- WEIGHTS
          if (any(M)) {
            W <- W + weps
            W <- W[!M] / sum(W[!M])
          }

          if (sum(!M) > 1) {
            SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples,
                               replace = TRUE, prob = W) 
          }
          else {
            SAMPLES <- rep( (1:nForecasts)[!M], nSamples)
          }

          tab <- rep(0, nForecasts)
          names(tab) <- members
          
          for (j in seq(along = tab)) tab[j] <- sum(SAMPLES == j)

          SAMPLES[] <- NA

          jj <- 0
          for (j in seq(along = tab)){
            nsamp <- tab[j]
            if (!nsamp) next
            SAMPLES[jj + 1:nsamp] <- rnorm( nsamp, MEAN[j], SD[j])
            jj <- jj + nsamp
          }

#        sampleMean[i] <- mean(SAMPLES) 
         sampleMedian[i] <- median(SAMPLES) 
  
# crps2 approximates a term that is quadratic in the number of members
         crps1  <- mean(abs(SAMPLES - obs[i])) 
         crps2 <-  mean(abs(diff(sample(SAMPLES))))
         crpsSim[i]  <- crps1 - crps2/2
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

