"ensembleForecasts.ensembleData" <-
function (x) 
{ 
 k <- attr(x, "ensembleSize")
 attr(x, "ensembleSize") <- NULL
 class(x) <- "data.frame"
 as.matrix(x[, 1:k])
}

