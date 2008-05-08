`julTOymdh` <-
function (julianDates, origin. = NULL, dropHour = NULL) 
{

 if (!is.null(orig <- attr( julianDates, "origin"))) {
   if (!is.null(origin.)) stop("origin is not uniquely specified")
   origin. <- orig
 }
 else if (is.null(origin.)) stop("origin is not specified")

 hour <- round(24*as.vector(julianDates - floor(julianDates)))
 x <- month.day.year( as.vector(round(julianDates)), origin. = origin.)

 if (is.null(dropHour)) {
   l <- attr(julianDates, "nchar")
   dropHour <- is.null(l) || l == 8
 }
 
 if (any(hour) || dropHour == FALSE) {
   x <- lapply(c(x[c("year", "month", "day")], list(hour = hour)), 
                 as.character)
 }
 else {
   x <- lapply(x[c("year", "month", "day")], as.character)
 }

 pad0mdh <- function(x) {
 pad0 <- function(x) if (nchar(x) == 2) x else paste("0", x, sep = "")
         as.vector(sapply(x, pad0))
         }

 x[-1] <- lapply(x[-1], pad0mdh)
 as.vector(apply(data.frame(x), 1, paste, collapse = ""))
}

