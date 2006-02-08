"ensembleObsNames.ensembleData" <-
function (x) 
{
 class(x) <- "data.frame" 
 names(x$obs)
}

