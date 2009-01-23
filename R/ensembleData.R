`ensembleData` <-
function(forecasts, dates = NULL, observations = NULL, ...,
         forecastHour, initializationTime, exchangeable = NULL)
{
 object <- match.call()
 object[[1]] <- as.name("cbind")
 object$forecastHour <- object$initializationTime <- object$exchangeable <- NULL
 object <- object[!unlist(lapply(object, is.null))]
 nam <- names(object)
 nam[2] <- nam[1]
 names(object) <- nam
 object <- eval(object, parent.frame())

 attr(object, "ensembleSize") <- if (is.null(dim(forecasts))) 1 else ncol(forecasts)
 if (missing(dates))
   stop("dates is missing")
 if (length(dates) == 1) {
   nObs <- if (is.null(dim(forecasts))) length(forecasts) else nrow(forecasts)
   dates <- as.factor(rep( as.character(dates), nObs))
 }
 if (missing(forecastHour)) {
   warning("data cannot be used for modeling unless forecast hour is specified")
  }
 else attr(object, "forecastHour") <- forecastHour
 if (missing(initializationTime)) {
   warning("initialization time not specified")
 }
 else  attr(object, "initializationTime") <- initializationTime
 attr(object, "exchangeable") <- exchangeable
 class(object) <- c("ensembleData", "data.frame")
 object
}

