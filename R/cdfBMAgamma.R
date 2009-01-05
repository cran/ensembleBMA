`cdfBMAgamma` <-
function (x, WEIGHTS, MEAN, VAR, offset = 0)
{
 RATE <- MEAN/VAR
 sum(WEIGHTS*pgamma(x,shape=RATE*MEAN,rate=RATE))-offset
}

