"ensembleForecasts.ensembleData" <-
function (x) 
{ 
 class(x) <- "data.frame"
 i <- charmatch("obs", dimnames(x)[[2]]) - 1
 as.matrix(x[, 1:i])
}

