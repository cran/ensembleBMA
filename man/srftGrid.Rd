\name{srftGrid}
\docType{data}
\alias{srftGrid}
\title{Gridded Surface Temperature Ensemble Forecasts}
\description{
This data set gives 48-hour forecasts of 2-m surface temperature and the
associated observations for the US Pacific Northwest initialized on 
February 14, 2004 OOZ and valid on Februrary 16, 2004 OOZ for a 102 by 99
grid. The ensemble forecasts come from an eight member version of the 
University of Washington Mesoscale Ensemble (UWME) that was operational 
in 2004. \cr
Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 10098 rows and 10 columns: \cr
CMCG        - forecasts from the CMCG member of the ensemble (numeric). \cr
ETA         - forecasts from the ETA member of the ensemble (numeric). \cr
GASP        - forecasts from the GASP member of the ensemble (numeric). \cr
GFS         - forecasts from the GFS member of the ensemble (numeric). \cr
CENT        - forecasts from the CENT member of the ensemble (numeric). \cr
JMA         - forecasts from the JMA member of the ensemble (numeric). \cr
NGAPS       - forecasts from the NGAPS member of the ensemble (numeric). \cr
TCWB        - forecasts from the TCWB member of the ensemble (numeric). \cr
UKMO        - forecasts from the UKMO member of the ensemble (numeric). \cr
latitude    - the latitude of each forecast/observation (numeric). \cr
longitude   - the longitude of each forecast/observation (numeric). \cr
}
\details{
This dataset is a subset of the data used in Berrocal et al. (2006),
provided for the purposes of testing. Typically forecasting would be
performed on much larger datasets.
}
\references{
V. J. Berrocal, A. E. Raftery and T. Gneiting,
Combining Spatial and Ensemble Information in Probabilistic Weather Forecasts,
Technical Report no. 495, Department of Statistics, University of
Washington, February 2006.

F. A. Eckel and C. F. Mass, 
Effective mesoscle, short-range ensemble forecasting,
\emph{Weather and Forecasting 20:328--350, 2005}.

E. P. Grimit and C. F. Mass, 
Initial results of a mesoscle short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.
}

\keyword{datasets}
