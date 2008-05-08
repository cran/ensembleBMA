`quantBMAnormal` <-
function(alpha, WEIGHTS, MEAN, SD)
{
  lower <- 0
  upper <- max(MEAN+6*SD)

  if (cdfBMAnormal(lower, WEIGHTS, MEAN, SD, 0) > alpha) return(NA)
  if (cdfBMAnormal(upper, WEIGHTS, MEAN, SD, 0) < alpha) return(NA)

  z <- uniroot(cdfBMAnormal, lower = lower, upper = upper,
               WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)

  z$root
}

