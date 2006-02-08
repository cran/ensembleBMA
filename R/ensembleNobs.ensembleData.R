"ensembleNobs.ensembleData" <-
function (x) 
{
 class(x) <- "data.frame" 
 length(x$obs)
}

