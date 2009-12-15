plot.fitBMAgamma <-
function(x, ensembleData, dates=NULL, ...) 
{

 exchangeable <- x$excchangeable
 
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 M <- matchEnsembleMembers(x,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 obs <- ensembleVerifObs(ensembleData)

 nForecasts <- ensembleSize(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- x$weights

 if (!all(Wmiss <- is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (x$varCoefs[1] + x$varCoefs[2]*f)^2
        
       fTrans <- sapply(f, powfun, power = x$power)

       MEAN <- apply(rbind(1, fTrans) * x$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       plotBMAgamma( WEIGHTS = W, MEAN = MEAN[!M], VAR = VAR[!M],
                     obs = obs[i], exchangeable = exchangeable)

    }

 }

 invisible()
}

