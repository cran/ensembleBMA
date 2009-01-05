`cdfBMAgamma0` <-
function (x, WEIGHTS, MEAN, VAR, PROB0, offset = 0)
{
 RATE <- MEAN/VAR
 sum(WEIGHTS*(PROB0+(1-PROB0)*pgamma(x,shape=MEAN*RATE,rate=RATE)))-offset
}

