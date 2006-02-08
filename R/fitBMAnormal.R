"fitBMAnormal" <-
function(ensembleData, control = controlBMAnormal())
{
  nObs <- ensembleNobs(ensembleData)
  ensMemNames <- ensembleMemberNames(ensembleData)
  nForecasts <- length(ensMemNames)

  if(is.null(sd <- control$start$sd)) sd <- sd(ensembleVerifObs(ensembleData))

  weights <- if (is.null(control$start$weights)) 1 else control$start$weights
  if (length(weights) == 1) weights <- rep(weights,nForecasts)
  weights <- pmax(weights,1.e-4)
  weights <- weights/sum(weights)
  if (!is.null(names(weights))) weights <- weights[ensMemNames]

  obs <- ensembleVerifObs(ensembleData)

# bias correction (constraints)

  switch(control$biasCorrection,
         regression = {
  meanFit <- apply(ensembleForecasts(ensembleData), 2, function(x, y) {
   components <- c("coefficients","fitted.values","residuals","model")
   lm(y~x)[components]},
    y = ensembleVerifObs(ensembleData))

  biasCoefs <- lapply(meanFit, function(x) x$coefficients)
  biasCoefs <- as.matrix(data.frame(biasCoefs))
  dimnames(biasCoefs) <- list(NULL, ensMemNames)
  bad <- biasCoefs[2,] < 0
  if (any(bad)) print("biasCoefs < 0") 

  MEAN <- lapply(meanFit, function(x) x$fitted.values)
  MEAN <- as.matrix(data.frame(MEAN))
  dimnames(MEAN) <- NULL

  RSQ <- lapply(meanFit, function(x) x$residuals)
  RSQ <- as.matrix(data.frame(RSQ))
  dimnames(RSQ) <- NULL
  RSQ <- RSQ^2
                       },
           additive = {
  intcpt <- apply(obs - ensembleForecasts(ensembleData), 
                   2, mean)
  biasCoefs <- rbind(intcpt,1)
  dimnames(biasCoefs) <- NULL
  MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 1, 
                FUN = "+", STATS = intcpt) 
  dimnames(MEAN) <- NULL
  RSQ <- ensembleVerifObs(ensembleData) - MEAN
  dimnames(RSQ) <- NULL
  RSQ <- RSQ^2
            },
           none = {
  biasCoefs <- matrix( c(0,1), 2, nForecasts) 
  MEAN <- ensembleForecasts(ensembleData)
  dimnames(MEAN) <- NULL
  RSQ <- obs - MEAN
  dimnames(RSQ) <- NULL
  RSQ <- RSQ^2
            }
 )

  z <- matrix(0, ncol=nForecasts, nrow=nObs)
  nIter <- 0
  loglik <- 0


  if (control$equalVariance)  sd <- rep(sd, nForecasts)

    while (TRUE) # EM
         {
          z <- sweep( dnorm(obs, mean=MEAN, sd=sd), 
                      MARGIN = 2, FUN = "*", STATS = weights)
#         z <- apply( z, 2, function(x) x /sum(x))
          zsum1 <- apply(z, 1, sum)
          z <- sweep( z, MARGIN = 1, FUN = "/", STATS = zsum1)

          old <- loglik
          loglik <- sum(log(zsum1))
 
          zsum2 <- apply(z, 2, sum)

          if (control$equalVariance) {
            sd <- sqrt(sum(z*RSQ)/sum(z))
           }
          else {
            sd <- sqrt(apply( z*RSQ, 2, sum)/zsum2)
          }
        
          weights <- zsum2/sum(zsum2)
        
          nIter <- nIter + 1    

          ERROR <- abs(loglik - old)/(1 + abs(loglik)) 
          if (nIter > 1 && ERROR < control$eps) break

          if (nIter >= control$maxIter) break
        }

 names(weights) <- ensMemNames

 structure(
  list(biasCoefs = biasCoefs, sd = sd, weights = weights, nIter = nIter),
  class = "fitBMAnormal")
}

