`ensembleData` <-
function(forecasts, caseLabels = NULL, memberLabels = NULL, 
         exchangeable = NULL,
         dates = NULL, observations = NULL, latitude = NULL,
         longitude = NULL, ...) 
{
 if (inherits(forecasts, "ensembleData")) class(forecasts) <- "data.frame"
 if (is.null(dim(forecasts))) forecasts <- as.matrix(forecasts)
 nObs <- nrow(forecasts)
 namesFor <- dimnames(forecasts)[[1]]
 if (is.null(memberLabels)) {
   memberLabels <- dimnames(forecasts)[[2]]
 }
 if (is.null(memberLabels)) {
   if (!is.null(exchangeable) && !is.null(names(exhangeable)))
     stop("exchangeable labels but no member labels")
   memberLabels <- as.character(1:ncol(forecasts))
   dimnames(forecasts) <- list(dimnames(forecasts)[[1]], memberLabels)
 }
 if (length(unique(memberLabels)) < length(memberLabels)) {
   stop("duplicated member labels")
 }
 if (!is.null(exchangeable)){
   if (length(exchangeable) != ncol(forecasts))
     stop("exchangeable specification incompatible with forecasts")
   namX <- names(exchangeable)
   if (!is.null(namX)) {
     if (length(unique(namX)) < length(namX)) {
       stop("duplicated exchangeable labels")
     }
     m <- match(namX, memberLabels, nomatch = 0)
     if (length(unique(m)) != length(m))
       stop("exchangeable incompatible with forecasts")
     exchangeable <- exchangeable[m]
   }
   else {
     names(exchangeable) <- memberLabels
   } 
 }
 forecasts <- as.data.frame(forecasts)
 ensembleSize <- ncol(forecasts)
 if (!is.null(dates)) {
    if (!all(dateCheck(dates))) stop("improperly specified date(s)")
    dates <- as.factor(dates)
  }
 object <- c(as.list(forecasts), list(observations = observations,
             dates = dates, 
             latitude = latitude, longitude = longitude), 
             list(...))
 object <- object[!sapply(object, is.null)]
 if (length(unique(sapply(object, length))) != 1) 
   stop("inputs unequal in length")
 namesObj <- names(object)
## names may change here e.g. avn/gfs to avn.gfs 
 object <- data.frame(object)
 namesObs <- names(observations)
 namesDates <- names(dates)
 namesFor <- dimnames(forecasts)[[1]]
 if (is.null(caseLabels)) {
    if (!is.null(namesFor)) {
      caseLabels <- namesFor
    }
    else if (!is.null(namesDates)) {
      caseLabels <- namesDates
    }
    else if (!is.null(namesObs)) {
      caseLabels <- namesObs
    }
 }
 dimnames(object) <- list(caseLabels, namesObj)

 attr(object, "ensembleSize") <- ensembleSize
 attr(object, "exchangeable") <- exchangeable
 class(object) <- c("ensembleData", "data.frame")
 object
}

