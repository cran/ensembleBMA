\name{slp}
\docType{data}
\alias{slp}
\title{48-hour sea-level pressure forecasts and observations for Puget Sound area}
\description{
This data set gives the 48-hour sea-level pressure forecasts and observations for the period from January 12, 2000 to June 30, 2000 at 0 hours GMT, using a 12-km grid resolution.
}
\format{
A matrix with nine columns and 496 rows. The columns are:

date - a numeric vector containing the date of each forecast/observation set, in the format YYYYMMDDHH.

lat - a numeric vector containing the latitude of each forecast/observation set.

lon - a numeric vector containing the latitude of each forecast/observation set.

F1 - a numeric vector containing the forecasts coming from the AVN member of the UW MM5.

F2 - a numeric vector containing the forecasts coming from the GEM member of the UW MM5.

F3 - a numeric vector containing the forecasts coming from the ETA member of the UW MM5.

F4 - a numeric vector containing the forecasts coming from the NGM member of the UW MM5.

F5 - a numeric vector containing the forecasts coming from the NOGAPS member of the UW MM5.

Y - a numeric vector containing the observed sea-level pressure.
}


\references{
Grimit, E.P. and C. F. Mass, "Initial results of a mesoscale short-range ensemble forecasting system over the Pacific Northwest." Weather and Forecasting, 2002, vol 17, pp 192-205
available at:
http://www.stat.washington.edu/MURI/papers/grimitmass2002.pdf

Raftery, A. E., T. Gneiting, F. Balabdaoui, & M. Polakowski, "Using Bayesian Model Averaging to calibrate forecast ensembles." Monthlly Weather Review, to appear, 2005.
earlier version available at:
http://www.stat.washington.edu/www/research/reports/2003/tr440.pdf
}
\keyword{datasets}
