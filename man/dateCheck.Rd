\name{dateCheck}
\alias{dateCheck}
\alias{getHH}
\title{
    Checks date format.
}
\description{
  Checks that the character form of a vector of dates conforms to
  YYYYMMDDHH or YYYYMMDD.
}
\usage{
  dateCheck(YYYYMMDDHH)
}
\arguments{
  \item{YYYYMMDDHH}{
    A character vector (or its factor equivalent) of dates which
    should be in the form YYYYMMDDHH or YYYYMMDD,
    in which YYYY specifies the year, MM the month,
    DD the day, and (optionally) HH the hour.
  }
}
\value{
  A logical vector indicating whether or not each element of YYYYMMDDHH
  has the correct format.
}
\details{
  If both YYYYMMDDHH and YYYYMMDD are present,
  the YYYYMMDD dates are assumed to be in error
  even if HH == 00 for all of the longer dates. \cr
  Requires the \code{chron} library.
}
\seealso{
  \code{\link{ymdhTOjul},\link{julTOymdh}}
}

\examples{
   dateCheck(c("2008043000", "20080431", "20080501"))
}
\keyword{chron}
% docclass is function
