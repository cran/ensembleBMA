`matchITandFH` <-
function(fit, ensembleData) {
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
UseMethod("matchITandFH")
}

`matchITandFH.default` <-
function( fit, ensembleData) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 naNULL <- function(x) {
            if (is.null(x)) return(-2^21)
            if (is.na(x)) return(-2^22)
            if (is.infinite(x)) return(-2^23)
            x
           }

 fitFH <- naNULL(attr(fit, "forecastHour"))
 fitIT <- naNULL(attr(fit, "initializationTime"))
 datFH <- naNULL(attr(fit, "forecastHour"))
 datIT <- naNULL(attr(fit, "initializationTime"))

 out <- TRUE

 if (fitFH != datFH & fitIT != datIT) {
   warning("forecast hour and initialization time inconsistent in data and fit\n")
   out <- FALSE
 }
 else if (fitFH != datFH) {
   warning("forecast hour inconsistent in data and fit\n")
   out <- FALSE
 }
 else if (fitIT != datIT) {
   warning("initialization time inconsistent in data and fit\n")
   out <- FALSE
 }

 out
}


