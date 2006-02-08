"gamma0BMAcdf" <-
function (x, offset, WEIGHTS, POP, MEAN, VAR)
{
  sum(WEIGHTS*(POP+(1-POP)*pgamma(rep(x,length(WEIGHTS)),
                                 shape=(MEAN^2/VAR),scale=VAR/MEAN)))-offset
}

