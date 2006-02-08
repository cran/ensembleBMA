"quantileForecastBMA" <-
function(fit, ensembleData, quantiles=0.5, dates=NULL, popData=NULL, ...) 
UseMethod("quantileForecastBMA")

