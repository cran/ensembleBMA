"gridForecastBMA.ensembleBMAgamma0" <-
function(object, gridData, quantiles = 0.5, date, popData = NULL, ...) 
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
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 k <- match(date, names(object$dateTable), nomatch=0)

 if (any(!k) || !length(k)) 
   stop("parameters not available for specified date")

 if (length(k) > 1) 
   stop("more than one date specified")

 ndim <- length(dimGrid <- dim(gridData)) 
 if (ndim!= 2 && ndim != 3) stop("gridData must be 2 or 3 dimensional")
 gridNames <- dimnames(gridData)

 if (ndim == 3) gridData <- apply( gridData, 3, as.vector)

 nGridPoints <- nrow(gridData)
 G <- matrix( NA, nGridPoints, length(quantiles))

 WEIGHTS <- object$weights[,k]
 
 for (i in 1:nGridPoints) {
    
       f <- gridData[i,]

       VAR <- object$varCoefs[1] + object$varCoefs[2]*f

       fTrans <- sapply(f, object$transformation)

       if (is.null(popData)) {
         POP <- sapply(apply(rbind( 1, fTrans, f==0) * object$popCoefs[,,k],
                                2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         POP <- sapply(apply(rbind( 1, fTrans, popi) * object$popCoefs[,,k],
                                2,sum), inverseLogit)
       }

       MEAN <- apply(rbind(1, fTrans)*object$biasCoefs[,,k], 2, sum)

       G[i,] <- sapply(sapply(quantiles,gamma0BMAquant,
                                WEIGHTS=WEIGHTS, POP=POP, MEAN=MEAN, VAR=VAR),
                         object$inverseTransformation)
 }

 if (ndim == 3) {
   G <- array( G, dimGrid)
   dimnames(G) <- list(gridNames[[1]], gridNames[[2]], as.character(quantiles))
 }
 else {
   dimnames(G) <- list(gridNames[[1]], as.character(quantiles))
 }

 G
}

