fitBMAgamma <-
function (ensembleData, control = controlBMAgamma(), exchangeable = NULL) 
{
    ZERO <- 1e-100
    powfun <- function(x, power) x^power
    if (is.null(exchangeable)) 
        exchangeable <- ensembleGroups(ensembleData)
    if (length(unique(exchangeable)) == length(exchangeable)) 
        exchangeable <- NULL
    if (!(nullEX <- is.null(exchangeable))) {
        namX <- as.character(exchangeable)
        uniqueX <- unique(namX)
        nX <- length(uniqueX)
    }
    maxIter <- control$maxIter
    tol <- eps <- control$tol
    ensembleData <- ensembleData[!dataNA(ensembleData, dates = FALSE), 
        ]
    nObs <- dataNobs(ensembleData)
    if (!nObs) 
        stop("no observations")
    obs <- dataVerifObs(ensembleData)
    if (is.null(startup <- dataStartupSpeed(ensembleData))) {
        if (is.null(control$startupSpeed)) 
            stop("default anemometer startup speed not specified")
        startup <- control$startupSpeed
    }
    if (length(startup) != nrow(ensembleData)) {
        startup <- rep(startup, length = nrow(ensembleData))
    }
    else if (length(startup) != 1) 
        stop("startup speed improperly specified")
    if (any(is.na(startup))) {
        if (is.null(control$startupSpeed)) 
            stop("default anemometer startup speed not specified")
        startup[is.na(startup)] <- control$startupSpeed
    }
    ensMemNames <- ensembleMembers(ensembleData)
    nForecasts <- length(ensMemNames)
    Y0 <- obs == 0
    if (sum(!Y0) < 2) 
        stop("less than 2 nonzero obs")
    ensembleData <- as.matrix(ensembleForecasts(ensembleData))
    ensembleData <- as.matrix(apply(ensembleData, 2, powfun, 
        power = control$power))
    obs <- as.vector(sapply(obs, powfun, power = control$power))
    lmFunc <- function(x, y) {
        beta0 <- min(y)
        x <- as.matrix(x)
        n <- ncol(x)
        x <- as.vector(x)
        nax <- is.na(x)
        x <- x[!nax]
        y <- rep(y, n)[!nax]
        if (all(!x)) {
            fit <- list(coefficients = c(mean(y), 0), fitted.values = rep(mean(y), 
                length(y)))
        }
        else {
            fit <- lm(y ~ x)
            coefs <- fit$coefficients
            if (coefs[1] <= 0) {
                coefs[1] <- beta0
                coefs[2] <- sum((y - beta0) * x)/sum(x * x)
                fit$coefficients <- coefs
                fit$fitted.values <- cbind(1, x) %*% coefs
            }
        }
        fit
    }
    meanFit <- lmFunc(ensembleData, obs)
    biasCoefs <- meanFit$coefficients
    meanVec <- as.vector(ensembleData)
    meanVec[!is.na(meanVec)] <- meanFit$fitted.values
    MEAN <- matrix(meanVec, nObs, nForecasts)
    miss <- is.na(ensembleData)
    Mzero <- miss[Y0, , drop = FALSE]
    Mnonz <- miss[!Y0, , drop = FALSE]
    completeDataLLmiss <- function(z, w, m, X, obs, startup) {
        objective <- function(par) {
            nObs <- length(obs)
            nFor <- ncol(X)
            miss <- is.na(X)
            Y0 <- obs == 0
            Mzero <- miss[Y0, , drop = FALSE]
            Mnonz <- miss[!Y0, , drop = FALSE]
            W <- matrix(w, nObs, nFor, byrow = TRUE)
            W[miss] <- 0
            W <- sweep(W, MARGIN = 1, FUN = "/", STATS = apply(W, 
                1, sum))
            v <- (par[1]^2 + (par[2]^2) * X)^2
            r <- m/v
            q <- array(NA, dim(z))
            q[Y0, ][!Mzero] <- log(pgamma(startup[Y0], shape = (r * 
                m)[Y0, , drop = FALSE][!Mzero], rate = r[Y0, 
                , drop = FALSE][!Mzero]))
            q[!Y0, ][!Mnonz] <- dgamma(matrix(obs, nObs, nFor)[!Y0, 
                , drop = FALSE][!Mnonz], shape = (r * m)[!Y0, 
                ][!Mnonz], rate = r[!Y0, , drop = FALSE][!Mnonz], 
                log = TRUE)
            include <- !miss & !(W == 0) & !(r == 0)
            -sum(z[include] * (q[include] + log(W[include])))
        }
        objective
    }
    varCoefs <- if (is.null(control$init$varCoefs)) 
        c(1, 1)
    else control$init$varCoefs
    varCoefs <- pmax(varCoefs, 1e-04)
    names(varCoefs) <- names(weights) <- NULL
    weights <- if (is.null(control$init$weights)) 
        1
    else control$init$weights
    if (length(weights) == 1) 
        weights <- rep(weights, nForecasts)
    weights <- weights/sum(weights)
    weights <- pmax(weights, 1e-04)
    weights <- weights/sum(weights)
    if (!is.null(names(weights))) 
        weights <- weights[ensMemNames]
    if (!nullEX) {
        for (labX in uniqueX) {
            I <- namX == labX
            weights[I] <- mean(weights[I])
        }
    }
    nIter <- 0
    z <- matrix(1/nForecasts, ncol = nForecasts, nrow = nObs)
    objold <- 0
    newLL <- 0
    while (TRUE) {
        VAR = (varCoefs[1] + varCoefs[2] * ensembleData)^2
        RATE <- MEAN/VAR
        SHAPE <- RATE * MEAN
        z[Y0, ][!Mzero] <- pgamma(startup[Y0], shape = SHAPE[Y0, 
            ][!Mzero], rate = RATE[Y0, ][!Mzero])
        z[!Y0, ][!Mnonz] <- dgamma(matrix(obs, nObs, nForecasts)[!Y0, 
            ][!Mnonz], shape = SHAPE[!Y0, ][!Mnonz], rate = RATE[!Y0, 
            ][!Mnonz], log = TRUE)
        zmax = apply(z[!Y0, , drop = F], 1, max, na.rm = TRUE)
        z[!Y0, ] <- exp(sweep(z[!Y0, , drop = F], MARGIN = 1, 
            FUN = "-", STATS = zmax))
        z <- sweep(z, MARGIN = 2, FUN = "*", STATS = weights)
        oldLL <- newLL
        newLL <- sum(zmax + log(apply(z[!Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        newLL <- newLL + sum(log(apply(z[Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        z <- z/apply(z, 1, sum, na.rm = TRUE)
        z[z < ZERO] <- 0
        wold <- weights
        zsum2 <- apply(z, 2, sum, na.rm = TRUE)
        weights <- zsum2/sum(zsum2)
        weights[weights < ZERO] <- 0
        if (!nullEX) {
            weights <- sapply(split(weights, namX), mean)[namX]
        }
        weps <- max(abs(wold - weights)/(1 + abs(weights)))
        fn <- completeDataLLmiss(z, weights, MEAN, ensembleData, 
            obs, startup)
        optimResult = if (is.null(control$optim.control)) {
            optim(sqrt(varCoefs), fn = fn, method = "BFGS")
        }
        else {
            optim(sqrt(varCoefs), fn = fn, method = "BFGS", control = control$optim.control)
        }
        if (optimResult$convergence) 
            warning("optim does not converge")
        varOld <- varCoefs
        varCoefs <- optimResult$par^2
        veps <- max(abs(varOld - varCoefs)/(1 + abs(varCoefs)))
        ERROR <- abs(objold - optimResult$value)/(1 + abs(optimResult$value))
        objold <- optimResult$value
        if (nIter > 0) {
            error <- abs(oldLL - newLL)/(1 + abs(newLL))
            if (error < eps) 
                break
        }
        nIter <- nIter + 1
        if (nIter >= maxIter) 
            break
    }
    if (nIter >= maxIter && error > eps) 
        warning("iteration limit reached")
    names(biasCoefs) <- NULL
    names(weights) <- ensMemNames
    startup <- unique(startup)
    if (length(startup) > 1) 
        startup <- NA
    structure(list(biasCoefs = biasCoefs, varCoefs = varCoefs, 
        weights = weights, nIter = nIter, loglikelihood = newLL, 
        power = control$power, startupSpeed = startup,
	call = match.call()), class = c("fitBMAgamma", "fitBMA"))
}
