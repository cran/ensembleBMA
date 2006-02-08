"quantileForecastBMA.fitBMAgamma0" <-
function(object, ensembleData, quantiles = 0.5, popData = NULL, ...) 
{

 inverseLogit <- function(x) {
# logit function safeguared against underflow and overflow
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  1/(1+exp(-x))
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  if (x >= log(.Machine$double.eps)) {
                    x <- exp(x)
                    x/(1+x)
                  }
                  else exp(x)
                }
                else 0
             }
            }

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (dim(popData) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)
 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(ensembleObsNames(ensembleData),as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

 for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       VAR <- object$varCoefs[1] + object$varCoefs[2]*f

       fTrans <- sapply(f, object$transformation)

       if (is.null(popData)) {
         POP <- sapply(apply(rbind( 1, fTrans, f==0) * object$popCoefs,
                                2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         POP <- sapply(apply(rbind( 1, fTrans, popi) * object$popCoefs,
                                2,sum), inverseLogit)
       }

       MEAN <- apply(rbind(1, fTrans)*object$biasCoefs, 2, sum)

       Q[i,] <- sapply(quantiles,gamma0BMAquant,WEIGHTS=object$weights,
                       POP=POP,MEAN=MEAN,VAR=VAR)
 }

  apply(Q, 2, object$inverseTransformation)
}

