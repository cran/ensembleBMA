\name{plotBMAforecast}
\alias{plotBMAforecast}
\alias{binGrid}
\title{
   Surface plots for BMA forecasts.
}
\description{
   Produces contour, image, or perspective plot of a BMA forecast
   using prediction on a grid.
}
\usage{
plotBMAforecast( forecast, longitude, latitude, nGrid = 65,
                 type = c("image", "contour", "persp"), \dots, 
                 interpolate = FALSE, span = 0.75, map = NULL) 
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
\item{map}{
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
  \code{\link{quantileForecastBMA}} 
}

\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
  data(srft)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")

  srftData <- ensembleData( forecasts = srft[,labels],
                            dates = srft$date, observations = srft$obs,
                            latitude = srft$lat, longitude = srft$lon)

\dontrun{
  bmaFit <- ensembleBMA( srftData, date = "2004012900",  
             model = "normal",trainingRule = list(length = 25, lag = 2))

  bmaForc <- quantileForecastBMA( bmaFit, srftData, date = "2004021400",  
                                  quantiles = c(.1, .5. .9))

  obs <- srftData$date == "2004012900"
  lat <- srftData$latitude[obs]
  lon <- srftData$longitude[obs]

  plotBMAforecast( bmaForc[,"0.5"], lat, lon, 
                   type = "contour", interpolate = TRUE)
  title("Median Forecast")

  plotBMAforecast( srftData$obs[obs], lat, lon, 
                  type = "contour", interpolate = TRUE)
  title("Observed Surface Temperature")

  data(srftGrid)

  memberLabels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
 
  srftGridData <- ensembleData(forecasts = srftGrid[,memberLabels],
      latitude = srftGrid[,"latitude"], longitude = srftGrid[,"longitude"])

  gridForc <- quantileForecastBMA( bmaFit, srftGridData, 
                    date = "2004021400", quantiles = c( .1, .5, .9))

  library(fields)

  plotBMAforecast(gridForc[,"0.5"],lon=srftGridData$lon,
     lat=srftGridData$lat,type="image",col=rev(rainbow(100,start=0,end=0.85)))
  title("Median Grid Forecast for Surface Temperature", cex = 0.5)

  probFreeze <- cdfBMA( bmaFit, srftGridData,  date = "2004021400", 
                             value = 273.15)

  plotBMAforecast(probFreeze, lon=srftGridData$lon, lat=srftGridData$lat,
                   type="image",col=gray((32:0)/32))
  title("Probability of Freezing", cex = 0.5)

}

}
\keyword{models}
% docclass is function
