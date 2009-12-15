ensembleValidDates.ensembleData <-
function (x) 
{ 
 class(x) <- "data.frame"
 as.character(x$dates)
}

