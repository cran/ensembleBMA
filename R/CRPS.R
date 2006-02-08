CRPS <- function(a,b,sigma,w,X,Y)
{
  #a couple of helper functions
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  absExp <- function(mu, sig)
  {
	sigma <- sig
	mean <- mu
	out <- (1/sqrt(pi))*sig*sqrt(2)*exp((-.5*mean^2)/(sig^2)) + mean*erf(.5*mean*sqrt(2)/sig)
	out
  }




  # Get the count of how many obs we have

  nmod <- ncol(X)

  # set up variance and sd vectors

  if(length(sigma)==1)
  {
    vars=rep(sigma^2, nmod)
  }
  else
  {
    vars=sigma^2
  }
  sds <- sqrt(vars)	

  # Create a matrix for the ensemble prediction

  means <- matrix(0,length(Y),nmod)
  maxFirst <- 0
  maxSecond <- 0
  i = 1
  while( i <= nmod )
  {
    means[,i] <- a[i] + b[i]*X[,i]
    i = i + 1
  }

  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))*sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.


  # Begin computing the first term in the CRPS formula.  This is a double sum since it is over w(i)*w(j) for all i and j.

  crps <- rep(0,length(Y))
  firSec <- matrix(0,length(Y),2)
  k = 1
  while( k <= nrow(means) )
  {
    firstSum <- 0
    secondSum <- 0
	
    # First get the first sum (the double sum of the CRPS)
	
    i = 1
    while( i <= nmod )
    {
      j = 1
      while( j <= nmod )
      {
        tvar <- vars[i] + vars[j]  # total variance
        tsd <- sqrt(tvar)          # total standard deviation
        tmean <- means[k,i] - means[k,j]
		
        firstSum <- firstSum + w[i]*w[j]*absExp(tmean,tsd)
        j = j + 1
      }
	
      i = i + 1
    }
	
    # Now get the second sum.  This one is only over all w(i) for all i.
	
    i = 1
    while( i <= nmod )
    {
      tvar <- vars[i]              # total variance
      tsd <- sqrt(tvar)            # total standard deviation
      tmean <- means[k,i] - Y[k]
      secondSum <- secondSum + w[i]*absExp(tmean,tsd)
      i = i + 1
    }
	
    # Using Szekely's expression for the CRPS, the first sum and second are put together
    # to compute the CRPS.

    crps[k]  <- -.5*firstSum + secondSum
    firSec[k,] <- c(.5*firstSum, secondSum)
		
    k = k + 1
  }

  out <- abs(mean(crps))
  out
}

