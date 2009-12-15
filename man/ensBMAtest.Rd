\name{ensBMAtest}
\docType{data}
\alias{ensBMAtest}
\alias{tempTestFit}
\alias{prcpTestFit}
\alias{winsTestFit}
\title{Ensemble BMA Test Data Set}
\description{
  This data set gives 48-hour forecasts for 2-m temperature,
  precipitation accumulated over the last 24 hours, and maximum wind speed
  at SeaTac (KSEA) and Portland (PDX) ariports in 2007/2008  initialized at 
  00 hours UTC using a 12km grid. The forecasts are based on an 8 member 
  version of the University of Washington mesoscale ensemble
  (Grimit and Mass 2002; Eckel and Mass 2005). \cr
}
\format{
A data frame with 66 rows and 34 columns: \cr
\code{idate} the initialization date of each forecast/observation, 
		format YYYYMMDDHH (categorical). \cr
\code{vdate} the validation date of each forecast/observation, 
		format YYYYMMDDHH (categorical). \cr
\code{latitude}  the latitude of each forecast/observation (numeric). \cr
\code{longitude}  the longitude of each forecast/observation (numeric). \cr
\code{longitude}  the elevation (in meters) above sea level (numeric). \cr
\code{station} weather station identifier (categorical). \cr
\code{network} weather network identifier (categorical). 
\code{*.gfs,*.cmcg,*.eta,*.gasp,*.jma,*.ngps,*.tcwb} 
forecasts from the 8 members of the ensemble (numeric).
\code{*.obs} observed values for the weather parameters.
The prefix \code{*} is one of \code{T2} for temperature,
    \code{PCP24} for precipitation, \code{MAXWSP10} for wind speed.\cr
}
\details{
  Temperature is given in Kelvin. \cr
  Precipitation amounts are quantized to hundredths of an inch. \cr
  Maximum wind speed is defined as the maximum of the hourly
  'instantaneous' wind speeds over the previous 18 hours, where an
  hourly 'instantaneous' wind speed is a 2-minute average from the
  period of two minutes before the hour to on the hour. \cr
  The wind speed observations are measured at 10-m above the ground and
  discretized when recorded by rounding to the
  nearest meter per second. \cr
This is a small dataset provided for the purposes of testing.
Typically forecasting would be performed on much larger datasets.
}
\references{
F. A. Eckel and C. F. Mass,
Effective mesoscale, short-range ensemble forecasting,
\emph{Weather and Forecasting 20:328--350, 2005}.

E. P. Grimit and C. F. Mass, 
Initial results of a mesoscale short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.
}
\examples{
\dontrun{# R check

  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

#----------------------------------------------------------------------------

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensMemNames, sep = ".")

  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)

  MAE( tempFit, tempTestData)
  CRPS( tempFit, tempTestData)

#----------------------------------------------------------------------------

  obs <- paste("PCP24","obs", sep = ".")
  ens <- paste("PCP24", ensMemNames, sep = ".")

  prcpTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs], 
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

  prcpTestFit <- ensembleBMAgamma0( prcpTestData, trainingDays = 30)

  MAE( prcpTestFit, prcpTestData)
  CRPS( prcpTestFit, prcpTestData)

#----------------------------------------------------------------------------

  obs <- paste("MAXWSP10","obs", sep = ".")
  ens <- paste("MAXWSP10", ensMemNames, sep = ".")

  winsTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

   winsTestFit <- ensembleBMAgamma(winsTestData, trainingDays = 30)

   MAE( winsTestFit, winsTestData)
   CRPS( winsTestFit, winsTestData)
}}
\keyword{datasets}
