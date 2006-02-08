"precipitationData" <-
function(forecasts, observations, dates, latitude = NULL, longitude = NULL,
         na.rm = TRUE, missingValues = NULL, labels = NULL)
{
 object <- ensembleData(forecasts, observations, dates, latitude, longitude,
                        na.rm, missingValues, labels)
 class(object) <- c("precipitationData", "ensembleData", "data.frame")
 object
}

