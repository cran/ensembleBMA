`ensembleBMA` <-
function(ensembleData, dates = NULL, trainingRule = list(length=NA, lag=NA), 
         control = NULL, model = NULL, warmStart = FALSE, minCRPS = NULL, 
         exchangeable = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 mc <- match.call()   
 mc$model <- NULL

 if (!is.null(model)) {
   MODELS <- c("normal","gamma0")
   m <- pmatch( model, MODELS, nomatch=0)
   model <- if (m) MODELS[m] else "?"
 }
 else stop("unspecified model")

  switch( model,
        "normal" = {
             mc[[1]] <- as.name("ensembleBMAnormal")
          },
        "gamma0" = {
             mc[[1]] <- as.name("ensembleBMAgamma0")
             if (!is.null(minCRPS)) {
               if (minCRPS) warning("minCRPS not available for gamma0 model")
               mc$minCRPS <- NULL
             }
         },
       stop("unrecognized model")
    )

 if (length(attr(ensembleData, "class")) > 2) {
   attr(ensembleData, "class") <- attr(ensembleData, "class")[-1]
   mc$ensembleData <- ensembleData
 }

 eval(mc, parent.frame())
}

