"fitBMA" <-
function(ensembleData, control = NULL, model = NULL, 
         exchangeable = NULL, popData = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
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
   if (!is.null(popData)) 
     warning("popData not used for normal models")
   mc$popData <- NULL
   mc[[1]] <- as.name("fitBMAnormal")
 }  
 else if (inherits(ensembleData, "precipitationData") || model == "gamma0") {
   mc[[1]] <- as.name("fitBMAgamma0")
 }
 else stop("unrecognized model")

 if (length(attr(ensembleData, "class")) > 2) {
   attr(ensembleData, "class") <- attr(ensembleData, "class")[-1]
   mc$ensembleData <- ensembleData
 }

 eval(mc, parent.frame())
}

