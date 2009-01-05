\name{ymdhTOjul}
\alias{ymdhTOjul}
\title{
    Convert to Julian dates.
}
\description{
  Converts YYYYMMDDHH or YYYYMMDD dates to Julian dates.
}
\usage{
 ymdhTOjul( YYYYMMDDHH, origin = c(month = 1, day = 1, year = 2000))
}
\arguments{
  \item{YYYYMMDDHH}{
    A character vector (or its factor equivalent) of dates in the form
    YYYYMMDDHH or YYYYMMDD, in which YYYY specifies the year, MM the
    month, DD the day, and (optionally) HH the hour.
  }
  \item{origin}{
    A named vector specifying the month, day, and year for the
    origin of the Julian dates. The default is
    \code{c(month = 1, day = 1, year = 2000)}.
  }
}
\value{
  A vector of Julian dates corresponding to YYYYMMDDHH.
  The vector has \code{"origin"} and \code{"dropHour"}
  attributes which give the origin for the Julian output and
  indicate whether or not the original format included the hour.
}
\details{
  Requires the \code{chron} library.
}
\seealso{
  \code{\link{dateCheck},\link{julTOymdh}}
}

\examples{
  data(ensBMAtest)

  julianVdates <- ymdhTOjul(ensBMAtest$vdate)
  all.equal( julTOymdh(julianVdates), as.character(ensBMAtest$vdate))

  all.equal( ymdhTOjul(ensBMAtest$idate), julianVdates-2)
}
\keyword{chron}
% docclass is function
