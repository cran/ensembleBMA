`fitBMAgamma0` <-
function(ensembleData, control = controlBMAgamma0(), exchangeable = NULL) 
{
  if (is.null(exchangeable)) exchangeable <- ensembleGroups(ensembleData)

  if (length(unique(exchangeable)) == length(exchangeable))
    exchangeable <- NULL

  if (!(nullX <- is.null(exchangeable))) {
    namX <- as.character(exchangeable)
    uniqueX <- unique(namX)
    nX <- length(uniqueX)
  }

  maxIter <- control$maxIter
  tol <- eps <- control$eps
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
  n0obs <- sum(Y0)

# untransformed weather data for variance model  

  ensembleData <- ensembleForecasts(ensembleData)

  Xvar <- ensembleData[!Y0, ]

##LM0 <- function(coefs, Xy, XX, R)(crossprod(R %*% coefs) - 2*sum(coefs*Xy))/2

##LM1 <- function(coefs, Xy, XX, R) XX %*% coefs - Xy

##LM2 <- function(coefs, Xy, XX, R) XX

  LM0 <- function(coefs, X, y, XX) {
    r <- y - X %*% coefs
    sum(r*r)/2
  }

  LM1 <- function(coefs, X, y, XX) {
    r <- y - X %*% coefs
    -crossprod(X,r)
  }

  LM2 <- function(coefs, X, y, XX) {
    XX
  }

  LMopt0 <-      function(z) {
## constrained opt for bias coefs with NA handling
# x is the cube root of the forecasts; y is the cube root of the obs
                           X <- model.matrix(z$formula, data = z$model)
                           coefs    <- z$coefficients 
                           coefs[1] <- max(coefs[1],0)
                           coefs[2] <- max(coefs[2],0)
                           opt <- nlminb(coefs, ob=LM0, gr=LM1, he=LM2,
                                         lower = c(0,0), upper =c(Inf,Inf),
                                     X = X, y = z$model$y,  XX = crossprod(X))
                           if (opt$convergence) print(opt$message)
                           coefs <- z$coefficients <- opt$par
                           z$fitted.values <- X %*% coefs
                           z
                          }

  LR0 <- function(coefs, X, y)
  {
    linearPredictor <- X %*% coefs
    z <- sapply(linearPredictor, function(x) {
                     if (x >= 0) {
                        if (-x >= log(.Machine$double.eps)) {
                          log(exp(-x)+1) + x 
                        }
                        else x
                      }
                     else {
                       if (x >= log(.Machine$double.eps)) {
                         log(1+exp(x))
                       }
                      else 0
                    }
             })
    sum(z - y*(linearPredictor))
  }

  LR1 <- function(coefs, X, y)
  {
    z <-  sapply(X %*% coefs,
                    function(x) {
                      if (x >= 0) {
                        if (-x >= log(.Machine$double.eps)) {
                          1/(1+exp(-x))
                        }
                        else 1
                      }
                     else {
                       if (x >= log(.Machine$double.xmin)) {
                         x <- exp(x)
                         x/(1+x)
                       }
                      else 0
                   }
             })
    -crossprod(X,y-z)
  }

  LR2 <- function(coefs, X, y)
  {
    d <- sapply(X %*% coefs,
                    function(x) {
                      if (x >= 0) {
                        if (-x >= log(.Machine$double.xmin)) {
                          x <- exp(-x)
                          x/(1+x)^2
                        }
                        else 0
                      }
                     else {
                       if (x >= log(.Machine$double.xmin)) {
                         x <- exp(x)
                         x/(1+x)^2
                       }
                      else 0
                     }
             })
    crossprod(sweep(X, MARGIN=1, FUN="*", STATS=sqrt(d)))
  }

  LRopt0 <-    function(z){
                       coefs <- z$coefficients
                       coefs[2] <- min(coefs[2],0)
                       coefs[3] <- max(coefs[3],0)
                       X <- model.matrix(z$formula, data = z$model)
                       opt <- nlminb(coefs, ob=LR0, gr=LR1, he=LR2, 
                                     lower = c(-Inf,-Inf,0), 
                                     upper =c(Inf,0,Inf),
                                     X = X, y = z$model$y)
                       if (opt$convergence)  print(opt$message)
                       coefs <- z$coefficients <- opt$par
                       z$fitted.values <- sapply(X %*% coefs, inverseLogit)
                       z
                      }
  
  ensembleData <- apply(ensembleData, 2, 
                       function(x) sapply(x,control$transformation))

  miss <- as.vector(as.matrix(is.na(ensembleData)))

  if (nullX) {
    prob0fit <- apply(ensembleData, 2, function(x, y0) {
    glm(y0~x+(x==0),family=binomial(logit), na.action = na.omit)}, y0 = Y0)

    prob0fit <- lapply( prob0fit, function(x, components) x[components],
       components = c("coefficients","fitted.values","model","formula"))

    temp <- lapply(prob0fit, function(x) x$fitted.values)

    prob0check <- unlist(lapply(prob0fit, function(z) {
                                coefs <- z$coefficients
                                coefs[2] <= 0 && coefs[3] >= 0}))

    if (any(!prob0check)) 
      prob0fit[!prob0check] <- lapply(prob0fit[!prob0check], LRopt0)

    prob0coefs <- lapply(prob0fit, function(x) x$coefficients)
    prob0coefs <- as.matrix(data.frame(prob0coefs))
    dimnames(prob0coefs) <- NULL

    PROB0 <- matrix(NA, nObs, nForecasts)
    PROB0[!miss] <- unlist(lapply(prob0fit, function(x) x$fitted.values))
  }
  else {

    logisticFunc <- function(x, y) {
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      y <- rep(y,n)
      glm(y ~ x + (x==0), na.action = na.omit)
    }

    prob0coefs <- matrix(NA, 3, nForecasts)
    dimnames(prob0coefs) <- NULL

    PROB0 <- matrix(NA, nObs, nForecasts)

    for (labX in uniqueX) {
       I <- namX == labX
       fit <- logisticFunc(ensembleData[, I, drop = F], Y0)
       coefs <- fit$coefficients
       if (any(!(fit$coefficients[2] <= 0 && fit$coefficients[3] >= 0))) {
         fit <- LRopt0(fit)
       }
       prob0coefs[,I] <- fit$coefficients
       miss <- is.na(ensembleData[, I, drop = F])
       miss <- as.vector(as.matrix(miss))
       PROB0[,I][!miss] <- fit$fitted.values
    }

  }

  PROB1 <- (1-PROB0)[!Y0,]
  PROB0 <- PROB0[Y0,]

  obs <- sapply(obs, function(x) sapply(x,control$transformation))
  obs <- obs[!Y0]
 
  nPrecip <- length(obs)

  ensembleData <- ensembleData[!Y0,]

  miss <- as.vector(as.matrix(is.na(ensembleData)))

# means determined as a bias-corrrrection step

  if (nullX) {
    meanFit <- apply(ensembleData, 2, function(x, y) {
              components <- c("coefficients","fitted.values","model","formula")
                     lm( y~x, na.action=na.omit)[components]}, y = obs)

    meanCheck1 <- unlist(lapply(meanFit, function(z) {
                                coefs <- z$coefficients
                                coefs[1] >= 0 && coefs[2] >= 0 
                               }))
    if (any(!meanCheck1)) 
      meanFit[!meanCheck1] <- lapply(meanFit[!meanCheck1], LMopt0)

     biasCoefs <- lapply(meanFit, function(x) x$coefficients)
     biasCoefs <- as.matrix(data.frame(biasCoefs))

     MEAN <- matrix(NA, nPrecip, nForecasts)
     MEAN[!miss] <- unlist(lapply(meanFit, function(x) x$fitted.values))
    }
  else {

    lmFunc <- function(x, y) {
      x <- as.matrix(x)
      n <- ncol(x)
      x <- as.vector(x)
      y <- rep(y,n)
      lm(y ~ x, na.action = na.omit)
    }

    biasCoefs <- matrix( NA, 2, nForecasts)

    MEAN <- matrix(NA, nPrecip, nForecasts)

    i <- 1
    for (labX in uniqueX) {
       I <- namX == labX
       fit <- lmFunc(ensembleData[, I, drop = F], obs)
       if (!(fit$coefficients[1] >=0 && fit$coefficients[2] >= 0)) 
         fit <- LMopt0(fit)
       biasCoefs[, I] <- fit$coefficients
       miss <- is.na(ensembleData[, I, drop = F])
       miss <- as.vector(as.matrix(miss))
       MEAN[,I][!miss] <- fit$fitted.values
    }

  }

  miss <- is.na(Xvar)

  gammaLoglikEM <- function(w, m, p1, X, Y)
{
  objective <- function(par)
 {
    v <- par[1]^2+(par[2]^2)*X

    g <- array(0,dim(v))
    rate <- m/v
    g <- dgamma(Y, shape=rate*m, rate=rate, log=TRUE)
    gmax <- max(g)
    g <- p1 * exp(g - gmax) # safeguard for over/underflow
    g  <- sweep(g, MARGIN=2, FUN= "*", STATS = w)
##  print(-sum(gmax+log(apply(g,1,sum))))
    -sum(gmax+log(apply(g,1,sum)))
  }
  objective
}

  gammaLoglikEMmiss <- function(w, m, p1, X, Y)
{
  objective <- function(par)
 {
    nObs <- length(Y)
    nFor <- ncol(X)

    G <- X
    miss <- is.na(X)

    W <- matrix( w, nObs, nFor, byrow = TRUE)
    W[miss] <- 0
    W <- sweep( W, MARGIN = 1, FUN = "/", STATS = apply(W, 1, sum))

    v <- par[1]^2+(par[2]^2)*X
    rate <- m/v
 
    G[!miss] <- dgamma(matrix(Y,nObs,nFor)[!miss],
                shape=(rate*m)[!miss], rate=rate[!miss], log=TRUE)

    gmax <- max(G, na.rm = TRUE)
    G <- p1 * exp(G - gmax) # safeguard for over/underflow
    G  <- G * W
    -sum(gmax+log(apply(G, 1, sum, na.rm = TRUE)), na.rm = TRUE)
  }
  objective
}

# gammaLoglikEMgrad <- function(par, w, m, p, X, Y)
 gammaLoglikEMgrad <- function(w, m, p1, X, Y)
{
  gradient <- function(par)
{
    v <- par[1]^2+(par[2]^2)*X

    g <- array(0,dim(v))
    alpha <- m^2/v
    g <- dgamma(Y, shape=alpha, rate=alpha/m, log=TRUE)
    gmax <- max(g)
    g <- p1*exp(g - gmax) # safeguard for over/underflow
    denom <- sum(gmax+log(apply(sweep(g,MARGIN=2,FUN= "*",STATS = w),1,sum)))
    g <- g * sweep(-digamma(alpha) + log(alpha-m) -
                   sweep(1/m, MARGIN=1, FUN="*", STATS = Y[!Y0]),
                   MARGIN=1, FUN="+", STATS = log(Y)+1) * alpha/v
    g  <- sweep(g, MARGIN=2, FUN= "*", STATS = w)
    -c(sum(apply(g,1,sum)),sum(apply(g*X,1,sum)))/denom
  }
  gradient
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
 
 for (i in 1:nEsteps) {

    z[!Y0,][!miss] <- dgamma(matrix(obs, nPrecip, nForecasts)[!miss], 
                             shape=SHAPE[!miss], rate=RATE[!miss], log=TRUE) 
    z[!Y0,] <- sweep(z[!Y0,], MARGIN=1, FUN="-", 
                     STATS=apply(z[!Y0,],1,max,na.rm=TRUE))  
    z[!Y0,] <- sweep(PROB1, MARGIN = 2, FUN="*", STATS=weights)*exp(z[!Y0,])
    z[Y0,] <- sweep(PROB0, MARGIN=2, FUN="*", STATS=weights)
 
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

} 

      fn <- gammaLoglikEMmiss(weights, MEAN, PROB1, Xvar, obs)
##     gr <- gammaLoglikEMgrad(weights, MEAN, PROB1, Xvar, obs)
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

    nIter <- nIter + 1

    if (nIter > 1 & ERROR < eps) break
    if (nIter >= maxIter) break
  }

## if (nIter >= maxIter && ERROR >= eps && max(c(veps,weps)) >= tol)
  if (nIter >= maxIter && ERROR >= eps)
    warning("iteration limit reached")

 dimnames(biasCoefs) <- list(NULL, ensMemNames)
 dimnames(prob0coefs) <- list(NULL, ensMemNames)
 names(weights) <- ensMemNames

  structure(
  list(prob0coefs = prob0coefs, biasCoefs = biasCoefs, varCoefs = varCoefs,
       weights = weights, nIter = nIter, 
       transformation = control$transformation,
       inverseTransformation = control$inverseTransformation),
       class = "fitBMAgamma0")
}

