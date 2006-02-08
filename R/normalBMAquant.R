"normalBMAquant" <-
function(alpha, WEIGHTS, MEAN, SD)
{
  lower <- 0
  upper <- max(MEAN+6*SD)

  if (normalBMAcdf(lower, 0, WEIGHTS, MEAN, SD) > alpha) return(NA)
  if (normalBMAcdf(upper, 0, WEIGHTS, MEAN, SD) < alpha) return(NA)

  z <- uniroot(normalBMAcdf, lower = lower, upper = upper,
               offset=alpha, WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD)

  z$root
}

