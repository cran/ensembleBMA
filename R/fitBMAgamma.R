`fitBMAgamma` <-
function(ensembleData, control = controlBMAgamma(), exchangeable = NULL) 
{

  powfun <- function(x,power) x^power

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
  nEsteps <- control$nEsteps

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

  if (sum(!Y0) < 2) stop("less than 2 nonzero obs")

# untransformed weather data for variance model  

  ensembleData <- ensembleForecasts(ensembleData)
  
  ensembleData <- as.matrix(apply(ensembleData, 2, 
                                  powfun, power = control$power))
  obs <- as.vector(sapply(obs, powfun, power = control$power))

# means determined as a bias-corrrrection step
# in this case, means are shared

  components <- c("coefficients","fitted.values","model")

  meanVec <- rep(is.na(obs),nForecasts) | as.vector(is.na(ensembleData))
  meanVec[as.logical(meanVec)] <- NA
  meanFit1 <- lm(rep(obs,nForecasts) ~ as.vector(ensembleData), 
           na.action = na.omit)[components]

  biasCoefs <- meanFit1$coef
  if (any(neg <- biasCoefs < 0)) {
    cat("bias coefs < 0")
    biasCoefs[neg] <- 0
  }

  meanVec[!is.na(meanVec)] <- meanFit1$fitted.values
 
  MEAN <- matrix(meanVec, nObs, nForecasts)

  miss <- is.na(ensembleData)

  Mzero <- miss[Y0,,drop=FALSE]
  Mnonz <- miss[!Y0,,drop=FALSE]

  gammaLoglikEM <- function(w, m, p1, X, Y)
{
  objective <- function(par)
 {
    v <- par[1]^2+(par[2]^2)*X

    g <- array(0,dim(v))
    rate <- m/v
    g <- dgamma(Y, shape=rate*m, rate=rate, log=TRUE)
    gmax <- apply(g, 1, max)
    g <- p1 * exp(g - gmax) # safeguard for over/underflow
    g  <- sweep(g, MARGIN=2, FUN= "*", STATS = w)
##  print(-sum(gmax+log(apply(g,1,sum))))
    -sum(gmax+log(apply(g,1,sum)))
  }
  objective
}

  gammaLoglikEMmiss <- function(w, m, X, Y)
{
  objective <- function(par)
 {
    nObs <- length(Y)
    nFor <- ncol(X)

    G <- X

    miss <- is.na(X)
    Y0 <- Y == 0

    Mzero <- miss[Y0,,drop=FALSE]
    Mnonz <- miss[!Y0,,drop=FALSE]

    W <- matrix( w, nObs, nFor, byrow = TRUE)
    W[miss] <- 0
    W <- sweep( W, MARGIN = 1, FUN = "/", STATS = apply(W, 1, sum))

    v <- (par[1]^2+(par[2]^2)*X)^2
    rate <- m/v
 
    G[Y0,][!Mzero] <- pgamma(1, shape=(rate*m)[Y0,][!Mzero], 
                                rate=rate[Y0,][!Mzero])
    G[!Y0][!Mnonz] <- dgamma(matrix(Y,nObs,nFor)[!Y0,][!Mnonz],
                shape=(rate*m)[!Y0,][!Mnonz], rate=rate[!Y0,][!Mnonz], log=TRUE)
    gmax <- rep(0,nrow(G))
    gmax[!Y0] <- apply( G[!Y0,], 1, max, na.rm = TRUE)
    G[!Y0,] <- exp(G[!Y0,] - gmax[!Y0]) # safeguard for over/underflow
    G  <- G * W
    -sum(gmax+log(apply(G, 1, sum, na.rm = TRUE)), na.rm = TRUE)
  }
  objective
}

  varCoefs <- if(is.null(control$start$varCoefs)) c(1,1) else control$start$varCoefs
  varCoefs <- pmax(varCoefs,1.e-4)

  names(varCoefs) <- names(weights) <- NULL

  # set all latent variables equal initially, 
  # and "new weights" (used later to compare changes in weights) to zero

  weights <- if (is.null(control$start$weights)) 1 else control$start$weights
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
    VAR= (varCoefs[1]+varCoefs[2]*ensembleData)^2

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

 for (i in 1:nEsteps) {

    z[Y0,][!Mzero] <- pgamma(1, shape=SHAPE[Y0,][!Mzero], rate=RATE[Y0,][!Mzero])

    z[!Y0,][!Mnonz] <- dgamma(matrix(obs, nObs, nForecasts)[!Y0,][!Mnonz], 

                  shape=SHAPE[!Y0,][!Mnonz], rate=RATE[!Y0,][!Mnonz], log=TRUE)

    zmax = apply( z[!Y0,], 1, max, na.rm=TRUE) 
    z[!Y0,] <- exp(sweep( z[!Y0,], MARGIN=1, FUN="-", STATS=zmax))

    z <- sweep( z, MARGIN = 2, FUN = "*", STATS = weights)

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

    if (nIter < 5) break
} 

      fn <- gammaLoglikEMmiss(weights, MEAN, ensembleData, obs)
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

 names(biasCoefs) <- NULL
 names(weights) <- ensMemNames

  structure(
  list(biasCoefs = biasCoefs, varCoefs = varCoefs,
       weights = weights, nIter = nIter, loglikelihood = newLL,
       power = control$power), class = "fitBMAgamma")
}

