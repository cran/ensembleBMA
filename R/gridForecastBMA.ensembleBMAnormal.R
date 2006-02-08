"gridForecastBMA.ensembleBMAnormal" <-
function(object, gridData, quantiles = 0.5, date, ...) 
{
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

 nForecasts <- dimGrid[ndim]

 SD <- if (!is.null(dim(object$sd))) {
         object$sd[,k] 
       }
       else {
         rep(object$sd[k], nForecasts)
       }
 
 WEIGHTS <- object$weights[,k]

 for (i in 1:nGridPoints) {
    
       f <- gridData[i,]

       MEAN <- apply(rbind(1, f)*object$biasCoefs[,,k], 2, sum)

       G[i, ] <- sapply(quantiles,normalBMAquant,
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

