"ensembleData" <-
function(forecasts, observations, dates, ..., na.action = "remove", 
          missingValues = NULL, labels = NULL) 
{
 auxList <- list(...)
 nObs <- length(observations)
 if (length(dates) != nObs) 
   stop("number of dates does not equal number of observations")
 dates <- as.factor(dates)
 if (is.null(dim(forecasts))) forecasts <- as.matrix(forecasts)
 if (nrow(forecasts) != nObs) 
   stop("number of observations does not equal number of forecasts")
 if (is.null(forcNames <- dimnames(forecasts)[[2]])) {
   forcNames <- paste("forc", 1:ncol(forecasts), sep = "")
 }
 object <- cbind(forecasts,observations,dates)
 NAMES <- c(forcNames, "observations", "dates")
#ord <- order(as.numeric(dates))
#object <- as.data.frame(object[ord,])
 namesObs <- names(observations)
 namesFor <- dimnames(forecasts)[[1]]
 if (is.null(labels)) {
    if (!is.null(namesObs <- names(observations))) {
      labels <- namesObs
    }
    else if (!is.null(namesFor <- dimnames(forecasts)[[1]])) {
      labels <- namesFor
    }
    else labels <- names(dates)
 }
 dimnames(object) <- list(labels, NAMES)
 object <- c(object, auxList)
 object <- as.data.frame(object)
 if (!is.null(missingValues)) {
     apply(object, 2, function(x, values) {
            x[match(values, x, nomatch=0)] <- NA
            }, values = as.vector(missingValues))
 }
 switch(na.action,
 remove= {
   object <- object[apply(object[,NAMES], 1, function(x) !any(is.na(x))),]
 },
 {})

 class(object) <- c("ensembleData", "data.frame")
 object
}

