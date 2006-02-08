"normalBMAquant" <-
function(alpha, WEIGHTS, MEAN, SD)
{
  lower <- 0
  upper <- max(MEAN+6*SD)

  if (normalBMAcdf(lower, WEIGHTS, MEAN, SD, 0) > alpha) return(NA)
  if (normalBMAcdf(upper, WEIGHTS, MEAN, SD, 0) < alpha) return(NA)

  z <- uniroot(normalBMAcdf, lower = lower, upper = upper,
               WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)

  z$root
}

