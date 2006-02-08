"fitBMAgamma0" <-
function(ensembleData, control = controlBMAgamma0(), popData = NULL) 
{
  maxIter <- control$maxIter
  tol <- eps <- control$eps
  nEsteps <- control$nEsteps

  ensMemNames <- ensembleMemberNames(ensembleData)
  nForecasts <- length(ensMemNames)

  obs <- ensembleVerifObs(ensembleData)
  nObs <- length(obs)
  Y0 <- obs == 0

# untransformed weather data for variance model  

  ensembleData <- ensembleForecasts(ensembleData)

  LM0 <- function(coefs, Xy, XX, R)(crossprod(R %*% coefs) - 2*sum(coefs*Xy))/2

  LM1 <- function(coefs, Xy, XX, R) XX %*% coefs - Xy

  LM2 <- function(coefs, Xy, XX, R) XX


 inverseLogit <- function(x) {
# logit function safeguared against underflow and overflow
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  x <- exp(-x)
                  1/(1+x)
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

    popFit <- apply(ensembleData, 2, function(x, y0) {
    # x is a transformation of the forecast
    glm(y0~x+(x==0),family=binomial(logit))}, y0 = Y0)

  if (is.null(popData)) {

    popFit <- lapply( popFit, function(x, components) x[components],
                        components = c("coefficients","fitted.values","model"))

    prob0check <- unlist(lapply(popFit, function(z) {
                                coefs <- z$coefficients
                                coefs[2] <= 0 && coefs[3] >= 0
                                }))
  if (any(!prob0check)) {
    popFit[!prob0check] <- lapply(popFit[!prob0check],
   function(z)  {
     coefs <- z$coefficients
     coefs[2] <- min(coefs[2],0)
     coefs[3] <- max(coefs[3],0)
     X <- model.matrix(y0~x+(x==0), data = z$model)
     opt <- nlminb(coefs, ob=LR0, gr=LR1, he=LR2, 
                   lower = c(-Inf,-Inf,0), upper =c(Inf,0,Inf),
                   X = X, y = z$model$y)
     if (opt$convergence)  print(opt$message)
     coefs <- z$coefficients <- opt$par
     z$fitted.values <- sapply(X %*% coefs, inverseLogit)
     z
          })
  }

  }
  else {

    if (!is.null(dim(popData))) {
      if (length(dim(popData)) == 2) {
        popData <- list(popData)
      }
     else {
       popData <- apply(popData, 3, list)
     }
    }
    
    for (j in 1:nForecasts) {
       popFit[[j]] <- update(popFit[[j]], . ~ . - (x == 0))
       popi <- do.call("cbind", lapply( popData, function(x,j) x[,j], j = j))
       popFit[[j]] <- update(popFit[[j]], . ~ . + popi)
    }

    popFit <- lapply( popFit, function(x, components) x[components],
                        components = c("coefficients","fitted.values","model"))

  }

  popCoefs <- lapply(popFit, function(x) x$coefficients)
  popCoefs <- as.matrix(data.frame(popCoefs))
  dimnames(popCoefs) <- NULL

  POP <- lapply(popFit, function(x) x$fitted.values)
  POP <- as.matrix(data.frame(POP))
  dimnames(POP) <- NULL
  PROB1 <- (1-POP)[!Y0,]
  POP <- POP[Y0,]

  Xvar <- ensembleData <- ensembleData[!Y0, ]

  obs <- sapply(obs, function(x) sapply(x,control$transformation))
  obs <- obs[!Y0]
  
  ensembleData <- apply(ensembleData, 2, 
                       function(x) sapply(x,control$transformation))


# means determined as a bias-corrrrection step
  meanFit <- apply(ensembleData, 2, function(x, y) {
   components <- c("coefficients","fitted.values","model")
   lm(y~x)[components]}, y = obs)

  meanCheck1 <- unlist(lapply(meanFit, function(z) {
                              coefs <- z$coefficients
                              coefs[1] >= 0 && coefs[2] >= 0 
                              }))
  if (any(!meanCheck1)) {
    meanFit[!meanCheck1] <- lapply(meanFit[!meanCheck1],
     function(z)  {
         y <- z$model$y
# x is the cube root of the forecasts; y is the cube root of the obs
         X <- model.matrix(y~x, data = z$model)
         coefs    <- z$coefficients 
         coefs[1] <- max(coefs[2],0)
         coefs[2] <- max(coefs[2],0)
     opt <- nlminb(coefs, ob=LM0, gr=LM1, he=LM2,
                   lower = c(0,0), upper =c(Inf,Inf),
                   Xy = crossprod(X,y), XX = crossprod(X), 
                   R = qr.R(qr(X)))
     if (opt$convergence)  print(opt$message)
     coefs <- z$coefficients <- opt$par
     z$fitted.values <- X %*% coefs
     z
          })
  }

# chk <- c(sum(as.numeric(!prob0check)), sum(as.numeric(!meanCheck1)))

# print(chk)

  biasCoefs <- lapply(meanFit, function(x) x$coefficients)
  biasCoefs <- as.matrix(data.frame(biasCoefs))
  dimnames(biasCoefs) <- NULL

  MEAN <- lapply(meanFit, function(x) x$fitted.values)
  MEAN <- as.matrix(data.frame(MEAN))
  dimnames(MEAN) <- NULL

#  gammaLoglikEM <- function(par, w, m, p, X, Y)
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
    -sum(gmax+log(apply(g,1,sum)))
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

  # set all latent variables equal initially, 
  # and "new weights" (used later to compare changes in weights) to zero

  weights <- if (is.null(control$start$weights)) 1 else control$start$weights
  if (length(weights) == 1) weights <- rep(weights,nForecasts) 
  weights <- pmax(weights,1.e-4)
  weights <- weights/sum(weights)
  if (!is.null(names(weights))) weights <- weights[ensMemNames]

  names(varCoefs) <- names(weights) <- NULL

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

# nEsteps <- max(nIter + 1, 100)
 for (i in 1:nEsteps) {

    z[!Y0,] <- dgamma(obs, shape=SHAPE, rate=RATE, log=TRUE) 
    z[!Y0,] <- sweep(z[!Y0,], MARGIN=1, FUN="-", STATS=apply(z[!Y0,],1,max))  
    z[!Y0,] <- sweep(PROB1, MARGIN=2, FUN="*", STATS=weights)*exp(z[!Y0,])
    z[Y0,] <- sweep(POP, MARGIN=2, FUN="*", STATS=weights)
 
    # normalize the latent variables
    z <- z/apply(z, 1, sum)

    # calculate new weights based on latent variables
    wold <- weights
    weights <- apply(z, 2, sum)/nObs
    weps <- max(abs(wold - weights)/(1+abs(weights)))

} 

      fn <- gammaLoglikEM(weights, MEAN, PROB1, Xvar, obs)
#     gr <- gammaLoglikEMgrad(weights, MEAN, PROB1, Xvar, obs)
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
 dimnames(popCoefs) <- list(NULL, ensMemNames)
 names(weights) <- ensMemNames

 structure(
  list(popCoefs = popCoefs, biasCoefs = biasCoefs, varCoefs = varCoefs,
       weights = weights, nIter = nIter, 
       transformation = control$transformation,
       inverseTransformation = control$inverseTransformation),
       class = "fitBMAgamma0")
}

