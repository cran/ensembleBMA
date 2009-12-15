fitBMAgamma0 <-
function(ensembleData, control = controlBMAgamma0(), exchangeable = NULL) 
{
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

  if (is.null(exchangeable)) exchangeable <- ensembleGroups(ensembleData)

  if (length(unique(exchangeable)) == length(exchangeable))
    exchangeable <- NULL

  if (!(nullX <- is.null(exchangeable))) {
    namX <- as.character(exchangeable)
    uniqueX <- unique(namX)
    nX <- length(uniqueX)
  }

  maxIter <- control$maxIter
  tol <- eps <- control$tol

# remove instances missing all forecasts or obs

  M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
  M <- M | is.na(ensembleVerifObs(ensembleData))
  ensembleData <- ensembleData[!M,]
 
  if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

# nObs <- length(obs)
  nObs <- ensembleNobs(ensembleData)

  ensMemNames <- ensembleMemberLabels(ensembleData)
  nForecasts <- length(ensMemNames)

  Y0 <- obs == 0
  n0obs <- sum(Y0)
  cat(" ", nObs - n0obs)

  if (sum(!Y0) < 2) stop("less than 2 nonzero obs")

# untransformed weather data for variance model  

  ensembleData <- ensembleForecasts(ensembleData)
  
  ensembleData <- apply(ensembleData, 2, 
                       function(x) sapply(x, powfun, power = control$power))

  miss <- as.vector(as.matrix(is.na(ensembleData)))

    logisticFunc <- function(x, y) {
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      y <- rep(y,n)
      nax <- is.na(x)
      x <- x[!nax]      
      y <- y[!nax]      
      if (!all(x==0) && !all(x!=0)) {
        fit <- glm(y ~ x + (x==0), family = binomial(logit))
        coefs <- fit$coefficients
      }
      else if (!all(x == 0)) {
        fit <- glm(y ~ x, family = binomial(logit))
        coefs <- c(fit$coefficients, 0)
      }
      else {
        fit <- glm(y ~ 1, family = binomial(logit))
        coefs <- c(fit$coefficients, 0, 0)
      }
      coefs[is.na(coefs)] <- 0
      if (!all(coefs[2:3] == 0) && coefs[2] > 0 && coefs[3] < 0) {
#       bad <<- cbind(y=y,x=x)
        stop("PROB0 coefficients: constraint violation")
      }
      else if (coefs[2] > 0) {
        fit <- glm(y ~ (x==0), family = binomial(logit))
        coefs <- c(fit$coefficients[1], 0, fit$coefficients[2])
        if (coefs[3] <= 0) stop("PROB0 constraint violation 3")
      }
      else if (coefs[3] < 0) {
        fit <- glm(y ~ x, family = binomial(logit))
        coefs <- c(fit$coefficients, 0)
        if (coefs[2] >= 0) stop("PROB0 constraint violation 2")
      }
     fit$coefficients <- coefs
     fit
    }

if (any(Y0)) {
  if (nullX) {
    prob0fit <- apply(ensembleData, 2, logisticFunc, y = Y0)

    prob0coefs <- lapply(prob0fit, function(x) x$coefficients)
    prob0coefs <- as.matrix(data.frame(prob0coefs))
    dimnames(prob0coefs) <- NULL

    PROB0 <- matrix(NA, nObs, nForecasts)
    PROB0[!miss] <- unlist(lapply(prob0fit, function(x) x$fitted.values))
  }
  else {

    prob0coefs <- matrix(NA, 3, nForecasts)
    dimnames(prob0coefs) <- NULL

    PROB0 <- matrix(NA, nObs, nForecasts)

    for (labX in uniqueX) {
       I <- namX == labX
       fit <- logisticFunc(ensembleData[, I, drop = FALSE], Y0)
       prob0coefs[,I] <- fit$coefficients
       miss <- is.na(ensembleData[, I, drop = FALSE])
       miss <- as.vector(as.matrix(miss))
       PROB0[,I][!miss] <- fit$fitted.values
    }

    p0 <- apply(PROB0, 1, max, na.rm = TRUE) == 1
    if (any(p0 & !Y0))  stop("PROB0 == 1 for a nonzero obs")

  }

  PROB1 <- (1-PROB0)[!Y0,]
  PROB0 <- PROB0[Y0,]

 }
 else {
      PROB1 <- matrix(1, nObs, nForecasts)
      prob0coefs <- matrix(0, 3, nForecasts)
      dimnames(prob0coefs) <- NULL
 }

  Xvar <- ensembleData[!Y0, ]

  obs <- sapply(obs, function(x) sapply( x, powfun, power = control$power))
  obs <- obs[!Y0]
 
  nPrecip <- length(obs)

  ensembleData <- ensembleData[!Y0,]

  miss <- as.vector(as.matrix(is.na(ensembleData)))

# means determined as a bias-corrrrection step

    lmFunc <- function(x, y) {
      beta0 <- min(y)
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      nax <- is.na(x)
      x <- x[!nax]
      y <- rep(y,n)[!nax]
      if (all(!x)) {
        fit <- list(coefficients = c(mean(y), 0),
                    fitted.values = rep(mean(y), length(y)))
      }
      else {
        fit <- lm(y ~ x)
        coefs <- fit$coefficients
        if (coefs[1] <= 0) {
          coefs[1] <- beta0
          coefs[2] <- sum((y-beta0)*x)/sum(x*x)
          fit$coefficients <- coefs
          fit$fitted.values <- cbind(1,x) %*% coefs
       }
     }
     fit
    }

  if (nullX) {

    meanFit <- apply(ensembleData, 2, lmFunc, y = obs)

     biasCoefs <- lapply(meanFit, function(x) x$coefficients)
     biasCoefs <- as.matrix(data.frame(biasCoefs))

     MEAN <- matrix(NA, nPrecip, nForecasts)
     MEAN[!miss] <- unlist(lapply(meanFit, function(x) x$fitted.values))
    }
  else {

    biasCoefs <- matrix( NA, 2, nForecasts)

    MEAN <- matrix(NA, nPrecip, nForecasts)

    i <- 1
    for (labX in uniqueX) {
       I <- namX == labX
       fit <- lmFunc(ensembleData[, I, drop = FALSE], obs)
       biasCoefs[, I] <- fit$coefficients
       miss <- is.na(ensembleData[, I, drop = FALSE])
       miss <- as.vector(as.matrix(miss))
       MEAN[,I][!miss] <- fit$fitted.values
    }

  }

  miss <- is.na(Xvar)

  completeDataLLmiss <- function(z, w, m, p0, p1, X, obs, Y0)
{
  objective <- function(par)
 {
# the obs contain no 0 values
    nObs <- length(obs)
    nFor <- ncol(X)

    miss <- is.na(X)

    v <- par[1]^2+(par[2]^2)*X
    r <- m/v

    W <- matrix( w, nObs, length(w), byrow = TRUE)
    W[miss] <- 0
    W <- sweep( W, MARGIN = 1, FUN = "/", STATS = apply(W, 1, sum))

    q <- array(NA, dim(z))

    q[!Y0,][!miss] <- dgamma(matrix(obs, nPrecip, nForecasts)[!miss], 
                             shape=(r*m)[!miss], rate=r[!miss], log=TRUE) 
#   q[!Y0,] <- log(sweep( p1, MARGIN = 2, FUN="*", STATS=w))+q[!Y0,]
#   q[Y0,] <- sweep( p0, MARGIN=2, FUN="*", STATS=w)

    Wzero <- (p1 == 0) | (W == 0)
    include <- !miss & !Wzero
    
    -sum(z[!Y0,][include]*(q[!Y0,][include]+log(p1[include]*W[include])))
  }
  objective
}

  varCoefs <- if(is.null(control$init$varCoefs)) c(1,1) else control$init$varCoefs
  varCoefs <- pmax(varCoefs,1.e-4)

  names(varCoefs) <- names(weights) <- NULL

  # set all latent variables equal initially, 
  # and "new weights" (used later to compare changes in weights) to zero

  weights <- if (is.null(control$init$weights)) 1 else control$init$weights
  if (length(weights) == 1) weights <- rep(weights,nForecasts) 
  weights <- weights/sum(weights)
  weights <- pmax(weights,1.e-4)
  weights <- weights/sum(weights)
  if (!is.null(names(weights))) weights <- weights[ensMemNames]

  if (!nullX) {
    for (labX in uniqueX) {
      I <- namX == labX
      weights[I] <- mean(weights[I])
    }
  }

  nIter <- 0
  z <- matrix( 1/nForecasts, ncol=nForecasts, nrow=nObs)
  objold <- 0

 # main EM algorithm

 newLL <- 0

  while(TRUE)
  {
    VAR= varCoefs[1]+varCoefs[2]*Xvar

    # set latent variables as weight times probability non-zero times gamma pdf
    # at that Y, # with the mean and variance parameters coming from each model
    # note that if Y equals zero, dgamma returns NaN since we only want our 
    # dgamma values for non-zero Y, and at zero Y we want weight times 
    # probability zero, we can now use the NaN entries as indicators of where 
    # to set latent variables to weight times probability zero

#   z=t(w*t((1-p)*dgamma((Y^(1/3)), shape=(m^2/v), rate=m/v)))
#   z[Y==0]=t(w*t(p))[Y==0]

# changed to compute at only non-zero values of observations (Chris F 8/06)

 RATE <- MEAN/VAR
 SHAPE <- RATE*MEAN

    z[!Y0,][!miss] <- dgamma(matrix(obs, nPrecip, nForecasts)[!miss], 
                             shape=SHAPE[!miss], rate=RATE[!miss], log=TRUE) 
    zmax = apply( z[!Y0,], 1, max, na.rm=TRUE) 
    z[!Y0,] <- sweep( z[!Y0,], MARGIN=1, FUN="-", STATS=zmax)  
    z[!Y0,] <- sweep( PROB1, MARGIN = 2, FUN="*", STATS=weights)*exp(z[!Y0,])
    z[Y0,] <- sweep( PROB0, MARGIN=2, FUN="*", STATS=weights)

    oldLL <- newLL

##  if (nIter > 0) {
##    opt <-  sum( zmax + log( apply( z[!Y0,], 1, sum, na.rm = TRUE))) 
##     print(c(opt, optimResult$value))
##   }

    newLL <-  sum( zmax + log( apply( z[!Y0,], 1, sum, na.rm = TRUE))) 
    newLL <- newLL + sum( log( apply( z[Y0,], 1, sum, na.rm=TRUE)))
 
    # normalize the latent variables
    z <- z/apply(z, 1, sum, na.rm = TRUE)

    # calculate new weights based on latent variables
    wold <- weights
    zsum2 <- apply(z, 2, sum, na.rm = TRUE)
    weights <- zsum2/sum(zsum2)

    if (!nullX) {

        weights <- sapply(split(weights,namX),mean)[namX]
##      for (labX in uniqueX) {
##        I <- namX == labX
##        weights[I] <- mean(weights[I])
##      }
    }

    weps <- max(abs(wold - weights)/(1+abs(weights)))

      fn <- completeDataLLmiss(z, weights, MEAN, PROB0, PROB1, Xvar, obs, Y0)
      optimResult = optim(sqrt(varCoefs), fn=fn, method = "BFGS") 

      if (optimResult$convergence) warning("optim does not converge")
      varOld <- varCoefs
      varCoefs <- optimResult$par^2
      veps <- max(abs(varOld - varCoefs)/(1+abs(varCoefs)))
      ERROR <- abs(objold - optimResult$value)/(1 + abs(optimResult$value))
      objold <- optimResult$value

# calculate change from last iteration to this one, 
# and then set weights to the new values

#   if (weps < tol && veps < tol) break

    if (nIter > 0) {
      error <- abs(oldLL - newLL)/(1 + abs(newLL))
     if (error < eps) break
    }

    nIter <- nIter + 1
    if (nIter >= maxIter) break
  }

## if (nIter >= maxIter && ERROR >= eps && max(c(veps,weps)) >= tol)
  if (nIter >= maxIter && error > eps)
    warning("iteration limit reached")

 dimnames(biasCoefs) <- list(NULL, ensMemNames)
 dimnames(prob0coefs) <- list(NULL, ensMemNames)
 names(weights) <- ensMemNames

  structure(
  list(prob0coefs = prob0coefs, biasCoefs = biasCoefs, varCoefs = varCoefs,
       weights = weights, nIter = nIter, loglikelihood = newLL,
       power = control$power), class = "fitBMAgamma0")
}

