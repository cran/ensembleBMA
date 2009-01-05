\name{julTOymdh}
\alias{julTOymdh}
\title{
    Convert Julian dates to character format.
}
\description{
  Converts Julian dates to YYYYMMDDHH or YYYYMMDD character format.
}
\usage{
 julTOymdh( julianDates, origin = NULL, dropHour = NULL)
}
\arguments{
  \item{julianDates}{
     A numeric vector specifying Julian dates.
  }
  \item{origin}{
    A named vector specifying the month, day, and year for the
    origin of the Julian dates. The default is
    \code{c(month = 1, day = 1, year = 2000)}.
    The default is \code{attr(julianDates, "origin")} if it exists.
  }
  \item{dropHour}{
    A logical value indicating whether of not the hour information
    should be drop from the specifiation of the dates if none of the
    Julian dates are fractional.
    The default is \code{attr(julianDates, "dropHour")} if it exists.
  }
}
\value{
    A character vector or numeric equivalent of dates in the form
    YYYYMMDDHH or YYYYMMDD, in which YYYY specifies the year, MM the
    month, DD the day, and (optionally) HH the hour
    corresponding to the Julian input.
}
\details{
  Requires the \code{chron} library.
}
\seealso{
  \code{\link{dateCheck},\link{ymdhTOjul}}
}

\examples{
  data(ensBMAtest)

  julianIdates <- ymdhTOjul(ensBMAtest$idate)

  all.equal( julTOymdh(julianIdates), as.character(ensBMAtest$idate))

  all.equal( ymdhTOjul(ensBMAtest$vdate), julianIdates+2)
}
\keyword{chron}
% docclass is function
