`CRPS.default` <-
function(fit, ensembleData, nSamples=NULL, seed=NULL, dates=NULL, ...) 
{
 mc <- match.call()   
 mc[[1]] <- as.name("crps")
 colMeans(eval(mc,parent.frame()))
}

