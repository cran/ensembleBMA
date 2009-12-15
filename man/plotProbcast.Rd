\name{plotProbcast}
\alias{plotProbcast}
\alias{binGrid}
\title{
   Surface plots for forecast information.
}
\description{
   Produces contour, image, or perspective plot of a forecast
   using loess prediction on a grid.
}
\usage{
plotProbcast( forecast, longitude, latitude, nGrid = 65,
                 type = c("image", "contour", "persp"), \dots, 
                 interpolate = FALSE, span = 0.75, maps = NULL) 
}
\arguments{
\item{forecast}{
   Numeric vector of forecasts.
 }
 \item{longitude}{
   Numeric vector giving the longitude of each forecast location.
 }
\item{latitude}{
   Numeric vector giving the latitude of each forecast location.
 }
 \item{nGrid}{
   Number of grid points for \code{loess} interpolation. 
  (Binning and interpolation are done on an \code{nGrid} by \code{nGrid} grid).
 }
 \item{type}{
   A character string indicating the desired plot type.
   Should be one of either \code{"contour"}, \code{"image"}, or
   \code{"persp"}.
 }
 \item{\dots}{
   Additional arguments to be passed to the plotting method.
 }
 \item{interpolate}{
   A logical variable indicating whether or not a \code{loess}
   fit should be used to interpolate the data to points on 
   a grid. The default is to determine grid values by binning,
   rather than interpolation.
 }
\item{span}{
   Smoothing parameter for \code{loess} (used only when
   \code{interpolate = TRUE}).
   The default value is \code{0.75}, which is 
   the default for \code{loess}.
 }
\item{maps}{
   A logical value indicating whether or not to include
   a map outline. The default is to include an outline
   if \code{type = "image"} and the \code{fields} library 
   is loaded.
 }
}
\value{
  An image, contour, or perspective plot of the forecast.
}
\details{
  If the \code{fields} library is loaded, a legend (and optionally
  a map outline) will be included in image plots.
}
\seealso{
  \code{\link{quantileForecast}} 
}

\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
}

\examples{
  data(srft)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")

  srftData <- ensembleData( forecasts = srft[,labels],
                            dates = srft$date, observations = srft$obs,
                            latitude = srft$lat, longitude = srft$lon,
                            forecastHour = 48, initializationTime = "00")

\dontrun{ # R check

  bmaFit <- ensembleBMA( srftData, date = "2004012900", trainingDays = 25,
                         model = "normal")

  bmaForc <- quantileForecast( bmaFit, srftData, date = "2004021400",  
                                  quantiles = c(.1, .5. .9))

  obs <- srftData$date == "2004012900"
  lat <- srftData$latitude[obs]
  lon <- srftData$longitude[obs]

  plotProbcast( bmaForc[,"0.5"], lat, lon, 
                   type = "contour", interpolate = TRUE)
  title("Median Forecast")

  plotProbcast( srftData$obs[obs], lat, lon, 
                  type = "contour", interpolate = TRUE)
  title("Observed Surface Temperature")

  data(srftGrid)

  memberLabels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
 
  srftGridData <- ensembleData(forecasts = srftGrid[,memberLabels],
      latitude = srftGrid[,"latitude"], longitude = srftGrid[,"longitude"],
                            forecastHour = 48, initializationTime = "00")

  gridForc <- quantileForecast( bmaFit, srftGridData, 
                    date = "2004021400", quantiles = c( .1, .5, .9))

  library(fields)

  plotProbcast(gridForc[,"0.5"],lon=srftGridData$lon,
     lat=srftGridData$lat,type="image",col=rev(rainbow(100,start=0,end=0.85)))
  title("Median Grid Forecast for Surface Temperature", cex = 0.5)

  probFreeze <- cdf( bmaFit, srftGridData,  date = "2004021400", 
                             value = 273.15)

  plotProbcast(probFreeze, lon=srftGridData$lon, lat=srftGridData$lat,
                   type="image",col=gray((32:0)/32))
  title("Probability of Freezing", cex = 0.5)

}}

\keyword{models}
% docclass is function
