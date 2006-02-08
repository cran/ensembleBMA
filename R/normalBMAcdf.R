"normalBMAcdf" <-
function (x, WEIGHTS, MEAN, SD, offset = 0)
{
  sum(WEIGHTS*pnorm(x, mean = MEAN, sd = SD)) - offset
}

