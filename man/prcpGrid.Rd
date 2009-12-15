\name{prcpGrid}
\alias{prcpGrid}
\title{
   Gridded Ensemble Forecasts of Precipitation
}

\description{
This data set gives 48-hour forecasts of 24 hour accumulated precipitation
on a grid of locations in the US Pacific Northwest 
initialized on  January 13, 2003 OOZ and valid on January 15, 2003 OOZ. 
The ensemble forecasts come from a nine member version of the 
University of Washington Mesoscale Ensemble 
(Grimit and Mass 2002; Eckel and Mass 2005).
Precipitation amounts are quantized to hundredths of an inch. 
%Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 8188 rows and 11 columns: \cr
\code{avn/gfs,cent,cmcg,eta,gasp,jma,ngps,tcwb,ukmo}
 forecasts from the 9 members of the ensemble (numeric). \cr
\code{latitude} the latitude of each forecast (numeric). \cr
\code{longitude} the longitude of each forecast (numeric). \cr
}

\references{
E. P. Grimit and C. F. Mass, 
Initial results of a mesoscale short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.

F. A. Eckel and C. F. Mass, 
Effective mesoscale, short-range ensemble forecasting,
\emph{Weather and Forecasting 20:328--350, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter, 
\code{ensembleBMA}: An \code{R} Package for Probabilistic Ensemble Forecasting 
using Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2009.
}

\examples{
\dontrun{ # R check

  data(prcpGrid) 

  prcpGridData <- ensembleData(forecasts = prcpGrid[,1:9], 
                               latitude = prcpGrid[,"latitude"],
                               longitude = prcpGrid[,"longitude"],
                               forecastHour = 48,
                               initilaizationTime = "00")

  data(prcpFit)

# median forecast for Jan 15, 2003 at the grid points

  quantileForecast( prcpFit, prcpGridData, date = "20030115")
  
}}
\keyword{data}

% docclass is function
