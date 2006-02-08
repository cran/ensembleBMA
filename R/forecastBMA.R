"forecastBMA" <-
function(date, ensembleData, trainingRule = list(length = 30, lag = 2), 
         control = NULL, quantiles = 0.5, model = NULL, 
         minCRPS = NULL, popData = NULL)
{
 if (!inherits(ensembleData, "ensembleData")) stop("not an ensembleData object")
 mc <- match.call()   
 mc$model <- NULL

 if (!is.null(model)) {
   MODELS <- c("normal","gamma0")
   m <- pmatch( model, MODELS, nomatch=0)
   model <- if (m) MODELS[m] else "?"
 }
 else model <- "?"

 if (inherits(ensembleData, "temperatureData") || 
     inherits(ensembleData, "pressureData") ||
     inherits(ensembleData, "precipitationData") || 
     model == "normal") {
   if (!is.null(popData)) {
     warning("popData not used for normal models")
     mc$popData <- NULL
   }
   if (is.null(minCRPS)) minCRPS <- TRUE
   mc[[1]] <- as.name("forecastBMAnormal")
 }  
 else if (inherits(ensembleData, "precipitationData") || model == "gamma0") {
   if (!is.null(minCRPS) && minCRPS) 
     warning("minCRPS not available for gamma0 models")
   mc$minCRPS <- NULL
   mc[[1]] <- as.name("forecastBMAgamma0")
 }
 else stop("unrecognized model")

 if (length(attr(ensembleData, "class")) > 2) {
   attr(ensembleData, "class") <- attr(ensembleData, "class")[-1]
   mc$ensembleData <- ensembleData
 }

 eval(mc, parent.frame())
}

