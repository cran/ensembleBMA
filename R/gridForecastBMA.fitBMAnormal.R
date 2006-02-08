"gridForecastBMA.fitBMAnormal" <-
function(object, gridData, quantiles = 0.5, ...) 
{
 ndim <- length(dimGrid <- dim(gridData)) 
 if (ndim!= 2 && ndim != 3) stop("gridData must be 2 or 3 dimensional")
 gridNames <- dimnames(gridData)

 if (ndim == 3) gridData <- apply( gridData, 3, as.vector)

 nGridPoints <- nrow(gridData)
 G <- matrix( NA, nGridPoints, length(quantiles))
 
 nForecasts <- dimGrid[ndim]

 SD <- if (!is.null(dim(object$sd))) {
         object$sd
       }
       else {
         rep(object$sd, nForecasts)
       }
 
 WEIGHTS <- object$weights

 for (i in 1:nGridPoints) {
    
       f <- gridData[i,]

       MEAN <- apply(rbind(1, f)*object$biasCoefs, 2, sum)

       G[i,] <- sapply(quantiles,normalBMAquant,
                         WEIGHTS=WEIGHTS,MEAN=MEAN,SD=SD)
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

