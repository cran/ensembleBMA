"bmaquant" <- function (a,b,sigma,w,alpha,X,niter=14)
{
  # Sept 12, 2003 (Adrian). Modified June 1, 2004
  # McLean 06/18/04 - changed arguments to bmacdf() to be consistent with new version of bmacdf()
  # Find the alpha quantile of the K-component BMA mixture 
  #  using the bisection method.
  # Inputs:
  #  a          vector of K intercepts in the regression bias correction
  #  b          vector of K slopes in the regression bias correction
  #  sigma      vector of K standard deviations from the BMA fit
  #             (a,b,sig are all outputs of EM.normals)
  #  alpha	required quantile
  #  X		vector of K forecasts
  #  niter	number of iterations in the bisection method
  #             Default 14, which gives accuracy on the order of
  #             (length of plausible interval)/16000

  out = NULL
  if( length(alpha) > 1 )
  {
	k = length(alpha)
	out = rep(0,k-1)
	for(i in 1:(k-1))
	{
		out[i] = bmaquant(a,b,sigma,w,alpha[i],X,niter)
	}
	alpha = alpha[k]
  }

  # Initialize: Find lower and upper bounds
  lower <- min (a+b*X-3*sigma)
  upper <- max (a+b*X+3*sigma)
  Flower <- bmacdf (a,b,sigma,w,X,lower)
  Fupper <- bmacdf (a,b,sigma,w,X,upper)
  if (Flower>alpha || Fupper<alpha) return(NA)

  # Bisection method
  for (iter in 1:niter)
  {
    mid <- (lower+upper)/2
    Fmid <- bmacdf (a,b,sigma,w,X,mid)
    if (Fmid>alpha) upper <- mid
    if (Fmid<alpha) lower <- mid
  }
  out = c(out, mid)
  out
}

