fitBMAnormal <-
function (ensembleData, control = controlBMAnormal(), exchangeable = NULL) 
{
    ZERO <- 1e-100
    if (is.null(exchangeable)) 
        exchangeable <- ensembleGroups(ensembleData)
    if (length(unique(exchangeable)) == length(exchangeable)) 
        exchangeable <- NULL
    if (!(nullEX <- is.null(exchangeable))) {
        namEX <- as.character(exchangeable)
        uniqueEX <- unique(namEX)
        nEX <- length(uniqueEX)
        splitEX <- split(seq(along = exchangeable), exchangeable)
        matEX <- sapply(uniqueEX, function(i, x) {
            x <- as.numeric(x == i)
            x/sum(x)
        }, x = exchangeable)
        dimnames(matEX) <- list(NULL, uniqueEX)
    }
    nObs <- dataNobs(ensembleData)
    if (!nObs) 
        stop("no observations")
    ensembleData <- ensembleData[!dataNA(ensembleData, dates = FALSE), 
        ]
    nObs <- dataNobs(ensembleData)
    if (!nObs) 
        stop("no observations")
    obs <- dataVerifObs(ensembleData)
    ensMemNames <- ensembleMembers(ensembleData)
    nForecasts <- length(ensMemNames)
    if (is.null(sd <- control$init$sd)) 
        sd <- sd(obs)
    weights <- if (is.null(control$init$weights) || any(is.na(control$init$weights))) 
        1
    else control$init$weights
    if (length(weights) == 1) 
        weights <- rep(weights, nForecasts)
    weights <- weights/sum(weights)
    weights <- pmax(weights, 1e-04)
    weights <- weights/sum(weights)
    if (!is.null(names(weights))) 
        weights <- weights[ensMemNames]
    MISSmat <- is.na(ensembleForecasts(ensembleData))
    MEAN <- RSQ <- matrix(NA, nObs, nForecasts)
    dimnames(MEAN) <- dimnames(RSQ) <- NULL
    switch(control$biasCorrection, regression = {
        if (nullEX) {
            meanFit <- apply(ensembleForecasts(ensembleData), 
                2, function(x, y) {
                  components <- c("coefficients", "fitted.values", 
                    "residuals", "model")
                  lm(y ~ x, na.action = na.omit)[components]
                }, y = obs)
            biasCoefs <- lapply(meanFit, function(x) x$coefficients)
            biasCoefs <- as.matrix(data.frame(biasCoefs))
            dimnames(biasCoefs) <- list(NULL, ensMemNames)
            FIT <- lapply(meanFit, function(x) x$fitted.values)
            RES <- lapply(meanFit, function(x) x$residuals)
            for (i in 1:nForecasts) {
                I <- !MISSmat[, i]
                MEAN[I, i] <- FIT[[i]]
                RSQ[I, i] <- RES[[i]]
            }
            MEAN <- as.matrix(data.frame(MEAN))
            dimnames(MEAN) <- NULL
            RSQ <- as.matrix(data.frame(RSQ))
            dimnames(RSQ) <- NULL
        } else {
            lmFunc <- function(x, y) {
                x <- as.matrix(x)
                n <- ncol(x)
                x <- as.vector(x)
                y <- rep(y, n)
                components <- c("coefficients", "fitted.values", 
                  "residuals", "model")
                lm(y ~ x, na.action = na.omit)[components]
            }
            biasCoefs <- matrix(NA, 2, nForecasts)
            dimnames(biasCoefs) <- list(NULL, ensMemNames)
            for (labX in uniqueEX) {
                I <- namEX == labX
                fit <- lmFunc(ensembleForecasts(ensembleData)[, 
                  I, drop = FALSE], obs)
                biasCoefs[, I] <- fit$coefficients
                M <- !MISSmat[, I, drop = FALSE]
                MEAN[, I][M] <- fit$fitted.values
                RSQ[, I][M] <- fit$residuals
            }
        }
        bad <- biasCoefs[2, ] < 0
        if (any(bad)) print("biasCoefs < 0")
    }, additive = {
        if (nullEX) {
            intcpt <- apply(obs - ensembleForecasts(ensembleData), 
                2, mean, na.rm = TRUE)
            biasCoefs <- rbind(intcpt, 1)
            dimnames(biasCoefs) <- NULL
            MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 2, 
                FUN = "+", STATS = intcpt)
            dimnames(MEAN) <- NULL
            RSQ <- obs - MEAN
            dimnames(RSQ) <- NULL
        } else {
            intFunc <- function(x, y) {
                x <- as.matrix(x)
                n <- ncol(x)
                x <- as.vector(x)
                y <- rep(y, n)
                mean(y - x, na.rm = TRUE)
            }
            intcpt <- rep(NA, nForecasts)
            for (labX in uniqueEX) {
                I <- namEX == labX
                intcpt[I] <- intFunc(ensembleForecasts(ensembleData)[, 
                  I, drop = FALSE], obs)
            }
            biasCoefs <- rbind(intcpt, 1)
            dimnames(biasCoefs) <- NULL
            MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 2, 
                FUN = "+", STATS = intcpt)
            dimnames(MEAN) <- NULL
            RSQ <- obs - MEAN
            dimnames(RSQ) <- NULL
        }
    }, none = {
        biasCoefs <- matrix(c(0, 1), 2, nForecasts)
        MEAN <- ensembleForecasts(ensembleData)
        dimnames(MEAN) <- NULL
        RSQ <- obs - MEAN
        dimnames(RSQ) <- NULL
    })
    RSQ <- RSQ^2
    z <- matrix(0, ncol = nForecasts, nrow = nObs)
    nIter <- 0
    loglik <- 0
    if (control$equalVariance) 
        sd <- rep(sd, nForecasts)
    while (TRUE) {
        z <- as.matrix(data.frame(lapply(1:nForecasts, function(i, 
            y, mu, sd) dnorm(y, mean = mu[, i], sd = sd[i]), 
            y = obs, mu = MEAN, sd = rep(sd, length = nForecasts))))
        z <- sweep(z, MARGIN = 2, FUN = "*", STATS = weights)
        dimnames(z) <- list(dimnames(z)[[1]], ensMemNames)
        zsum1 <- apply(z, 1, sum, na.rm = TRUE)
        if (any(zsum1 == 0)) 
            stop("parameter estimation fails")
        z <- sweep(z, MARGIN = 1, FUN = "/", STATS = zsum1)
        z[z < ZERO] <- 0
        old <- loglik
        loglik <- sum(log(zsum1))
        zsum2 <- apply(z, 2, sum, na.rm = TRUE)
        weights <- zsum2/sum(zsum2)
        weights[weights < ZERO] <- 0
        if (nullEX) {
            if (control$equalVariance) {
                sd <- sqrt(sum(z * RSQ, na.rm = TRUE)/sum(z, 
                  na.rm = TRUE))
            }
            else {
                sd <- sqrt(apply(z * RSQ, 2, sum, na.rm = TRUE)/zsum2)
            }
        }
        else {
            weights <- drop(crossprod(weights, matEX)[, namEX])
            if (control$equalVariance) {
                sd <- sqrt(sum((z * RSQ) %*% matEX, na.rm = TRUE)/sum(z %*% 
                  matEX, na.rm = TRUE))
                z <- sweep(z, MARGIN = 2, FUN = "/", STATS = nEX)
            }
            else {
                sd <- sqrt(apply(z * RSQ, 2, sum, na.rm = TRUE)/zsum2)
                sd <- drop(crossprod(sd, matEX))[namEX]
            }
        }
        nIter <- nIter + 1
        ERROR <- abs(loglik - old)/(1 + abs(loglik))
        if (nIter > 1 && ERROR < control$tol) 
            break
        if (nIter >= control$maxIter) 
            break
    }
    names(weights) <- ensMemNames
    structure(list(biasCoefs = biasCoefs, sd = sd, weights = weights, 
        exchangeable = exchangeable,
	nIter = nIter, loglikelihood = loglik, call = match.call()), 
        class = c("fitBMAnormal", "fitBMA"))
}
