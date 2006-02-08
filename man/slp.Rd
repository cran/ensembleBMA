\name{slp}
\docType{data}
\alias{slp}
\title{Sea-Level Pressure Ensemble Forecasts and Observations}
\description{
This data set gives 48-hour sea-level pressure forecasts for the
US Pacific Northwest from a 5 member ensemble at the University of Washington 
(UW MM5) and the associated observations for five different stations
for the period from April 16, 2000 to June 30, 2000 
at 0 hours GMT, using a 12-km grid resolution. \cr
Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 205 rows and 11 columns: \cr
AVN         - forecasts from the AVN member of the ensemble (numeric). \cr
GEN         - forecasts from the GEM member of the ensemble (numeric). \cr
ETA         - forecasts from the ETA member of the ensemble (numeric). \cr
NGM         - forecasts from the NGM member of the ensemble (numeric). \cr
NOGAPS      - forecasts from the NOGAPS member of the ensemble (numeric). \cr
observation - the observed sea-level pressure (numeric). \cr
date        - the date of each forecast/observation set, \cr
		in the format YYYYMMDDHH (categorical). \cr
latitude    - the latitude of each forecast/observation (numeric). \cr
longitude   - the longitude of each forecast/observation (numeric). \cr
station     - weather station identifier (categorical). \cr
type        - weather station type (categorical). \cr
}
\details{
This dataset is a small subset of the data used in Raftery et al. (2005),
provided for the purposes of testing. Typically forecasting would be
performed on much larger datasets.
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

E. P. Grimit and C. F. Mass, 
Initial results of a mesoscle short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.
}
\keyword{datasets}
