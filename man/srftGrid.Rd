\name{srftGrid}
\docType{data}
\alias{srftGrid}
\title{Gridded Surface Temperature Ensemble Forecasts}
\description{
This data set gives 48-hour forecasts of 2-m surface temperature
on a grid of locations in the US Pacific Northwest 
initialized on January 29, 2004 00UTC and valid on January 31, 2004 00UTC.
The ensemble forecasts come from an eight member version of the 
University of Washington Mesoscale Ensemble 
(Grimit and Mass 2002; Eckel and Mass 2005).
Temperatures are measured in kelvins. \cr
Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 10098 rows and 10 columns: \cr
\code{CMCG,ETA,GASP,GFS,JMA,NGAPS,TCWB,UKMO}
 forecasts from the 8 members of the ensemble (numeric). \cr
\code{latitude} the latitude of each forecast (numeric). \cr
\code{longitude} the longitude of each forecast (numeric). \cr
}


\references{
F. A. Eckel and C. F. Mass, 
Effective mesoscale, short-range ensemble forecasting,
\emph{Weather and Forecasting 20:328--350, 2005}.

E. P. Grimit and C. F. Mass, 
Initial results of a mesoscale short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.

V. J. Berrocal, A. E. Raftery and T. Gneiting,
Combining spatial and ensemble information in probabilistic weather forecasts,
\emph{Monthly Weather Review 133:1386--1402, 2007}.

V. J. Berrocal, A. E. Raftery, T. Gneiting and R. C. Steed,
Probabilistic Weather Forecasting for Winter Road Maintenance,
\emph{Journal of the American Statistical Association, 2010 (to appear)}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
}
\examples{
\dontrun{ # R check
  data(srft)
  data(srftGrid)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")

  srftData <- ensembleData( forecasts = srft[ ,labels],
                            dates = srft$date,
                            observations = srft$obs,
                            latitude = srft$lat,
                            longitude = srft$lon,
                            forecastHour = 48,
                            initializationTime = "00")

  srftFit <- ensembleBMAnormal( srftData, date = "2004013100",
                                trainingDays = 25)


  srftGridData <- ensembleData( forecasts = srftGrid[ ,labels],
                                latitude = srftGrid$lat,
                                longitude = srftGrid$lon,
                                forecastHour = 48,
                                initializationTime = "00")

   CRPS( srtGridData, srftFit)
}}
\keyword{datasets}
