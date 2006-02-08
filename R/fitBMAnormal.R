"fitBMAnormal" <-
function(ensembleData, control = controlBMAnormal(),
         exchangeable = NULL)
{
  if (is.null(exchangeable)) exchangeable <- ensembleGroups(ensembleData)

  if (length(unique(exchangeable)) == length(exchangeable))
    exchangeable <- NULL
 
  if (!(nullX <- is.null(exchangeable))) {
    namEX <- as.character(exchangeable)
    uniqueEX <- unique(namEX)
    nEX <- length(uniqueEX)
    splitEX <- split(seq(along = exchangeable), exchangeable)
    matEX <- sapply(uniqueEX, function(i,x) {
                       x <- as.numeric(x == i);
                       x/sum(x)}, x = exchangeable)
    dimnames(matEX) <- list(NULL, uniqueEX)
  }

  nObs <- ensembleNobs(ensembleData)
  ensMemNames <- ensembleMemberLabels(ensembleData)
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
  if (nullX) {
    meanFit <- apply(ensembleForecasts(ensembleData), 2, function(x, y) {
    components <- c("coefficients","fitted.values","residuals","model")
    lm(y~x)[components]},
    y = ensembleVerifObs(ensembleData))

    biasCoefs <- lapply(meanFit, function(x) x$coefficients)
    biasCoefs <- as.matrix(data.frame(biasCoefs))
    dimnames(biasCoefs) <- list(NULL, ensMemNames)

    MEAN <- lapply(meanFit, function(x) x$fitted.values)
    MEAN <- as.matrix(data.frame(MEAN))
    dimnames(MEAN) <- NULL

    RSQ <- lapply(meanFit, function(x) x$residuals)
    RSQ <- as.matrix(data.frame(RSQ))
    dimnames(RSQ) <- NULL
  }
  else {

    lmFunc <- function(x, y) {
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      y <- rep(y,n)
      components <- c("coefficients","fitted.values","residuals","model")
      lm(y ~ x)[components]
    }

    biasCoefs <- matrix(NA,2,nForecasts)
    dimnames(biasCoefs) <- list(NULL, ensMemNames)

    MEAN <- RSQ <- matrix(NA,nObs,nForecasts)
    dimnames(MEAN) <- dimnames(RSQ) <- NULL

    for (labX in uniqueEX) {
       I <- namEX == labX
       fit <- lmFunc(ensembleForecasts(ensembleData)[, I, drop = F],
                     ensembleVerifObs(ensembleData))
       biasCoefs[, I] <- fit$coefficients
       MEAN[,I] <- fit$fitted.values
       RSQ[,I] <- fit$residuals
    }

  }

  bad <- biasCoefs[2,] < 0
  if (any(bad)) print("biasCoefs < 0") 
                     },
          additive = {
  if (nullX) {
    intcpt <- apply(obs - ensembleForecasts(ensembleData), 2, mean)
   
    biasCoefs <- rbind(intcpt,1)
    dimnames(biasCoefs) <- NULL
   
    MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 2, 
                  FUN = "+", STATS = intcpt) 
    dimnames(MEAN) <- NULL
    RSQ <- ensembleVerifObs(ensembleData) - MEAN
    dimnames(RSQ) <- NULL
  }
  else {

    intFunc <- function(x, y, EX) {
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      y <- rep(y,n)
      uEX <- unique(EX)
      nEX <- length(uEX)
      intcpt <- matrix(NA,length(y),nE)
      for (j in 1:nEX) {
         intcpt[,j] <- y - apply(x[,EX == uEX[j],drop=FALSE],2,mean)
      }
      apply(intcpt, 2, mean)
    }

    intcpt <- intFunc( ensembleForecasts(ensembleData), 
                       ensembleVerifObs(ensembleData), exchangeable)
   
    biasCoefs <- rbind(intcpt,1)
    dimnames(biasCoefs) <- NULL
   
    MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 2, 
                  FUN = "+", STATS = intcpt) 
    dimnames(MEAN) <- NULL

    RSQ <- ensembleVerifObs(ensembleData) - MEAN
    dimnames(RSQ) <- NULL
  }
                     },
              none = {
  biasCoefs <- matrix( c(0,1), 2, nForecasts) 
  MEAN <- ensembleForecasts(ensembleData)
  dimnames(MEAN) <- NULL
  RSQ <- obs - MEAN
  dimnames(RSQ) <- NULL
                     }
  )

  RSQ <- RSQ^2

  z <- matrix(0, ncol=nForecasts, nrow=nObs)
  nIter <- 0
  loglik <- 0

  if (control$equalVariance)  sd <- rep(sd, nForecasts)

  xSDeq <- function( z, RSQ, splitEX){
      zr <- zz <- 0
      RSQ <- z*RSQ
      for (I in splitEX) {
         zr <- zr + sum(apply( RSQ[, I, drop = FALSE], 1, mean))
         zz <- zz + sum(apply(z[, I, drop = FALSE], 1, mean))
      }
      sqrt(zr/zz)
  }

##  cat("\n")
    while (TRUE) # EM
         {
          z <- as.matrix(data.frame(lapply(1:nForecasts, 
                    function(i, y, mu, sd) dnorm(y, mean=mu[,i], sd=sd[i]), 
                      y = obs, mu = MEAN, sd = rep(sd, length = nForecasts))))

          z <- sweep( z, MARGIN = 2, FUN = "*", STATS = weights)
          dimnames(z) <- list(dimnames(z)[[1]], ensMemNames)

          zsum1 <- apply(z, 1, sum)

          z <- sweep( z, MARGIN = 1, FUN = "/", STATS = zsum1)

          old <- loglik
          loglik <- sum(log(zsum1))
 
          zsum2 <- apply(z, 2, sum)
 
          weights <- zsum2/sum(zsum2)

          if (nullX) {
            if (control$equalVariance) {
              sd <- sqrt(sum(z*RSQ)/sum(z))
            }
            else {
              sd <- sqrt(apply(z*RSQ,2,sum)/zsum2)
            }
          }
          else {

##          weights <- sapply(split(weights,namEX), mean)[namEX]
            weights <- drop(crossprod(weights,matEX))[namEX]
            if (control$equalVariance) {
              sd <- sqrt(sum((z*RSQ)%*%matEX)/sum(z%*%matEX))
            }
            else {
              sd[] <- sqrt(apply(z*RSQ,2,sum)/zsum2)
              sd <- drop(crossprod(sd,matEX))[namEX]
            }
          }
        
          nIter <- nIter + 1    

          ERROR <- abs(loglik - old)/(1 + abs(loglik)) 
##          cat("", nIter)
##          cat(" ", c(loglik, ERROR, min(weights)), "\n")
          if (nIter > 1 && ERROR < control$eps) break

          if (nIter >= control$maxIter) break
        }

 names(weights) <- ensMemNames

 structure(
  list(biasCoefs = biasCoefs, sd = sd, weights = weights, 
       exhangeable = exchangeable, nIter = nIter),
  class = "fitBMAnormal")
}

