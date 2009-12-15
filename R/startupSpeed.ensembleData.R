startupSpeed.ensembleData <-
function (x) 
{
 if (is.null(startup <- attr( x, "startupSpeed"))) {
   class(x) <- "data.frame" 
   x$startup
 }
 else startup
}

