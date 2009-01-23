`verifRank` <-
function(forecasts, observations) 
{
 apply((forecasts > observations),1,sum)
}

