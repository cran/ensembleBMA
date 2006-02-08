"normalBMAcdf" <-
function (x, offset, WEIGHTS, MEAN, SD)
{
  sum(WEIGHTS*pnorm(x, mean = MEAN, sd = SD)) - offset
}

