"gamma0BMAcdf" <-
function (x, WEIGHTS, PROB0, MEAN, VAR, offset = 0)
{
 sum(WEIGHTS*(PROB0+(1-PROB0)*pgamma(x,shape=(MEAN^2/VAR),scale=VAR/MEAN)))-offset
}

