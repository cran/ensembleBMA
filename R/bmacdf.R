

"bmacdf" <- function (a,b,sigma,w,X,Y)
{
  # Sept 12, 2003 (Adrian)
  # McLean 06/28/04 - changed hard-coded 8 dimensions to arbitrary number of dimensions
  # Find the BMA cdf at x for the K-component BMA mixture 
  # Inputs:
  #  a          vector of K intercepts in the regression bias correction
  #  b          vector of K slopes in the regression bias correction
  #  sigma        vector of K standard deviations from the BMA fit
  #  w          vector of K weights from the BMA fit
  #  Y		value of temperature or SLP at which cdf is required
  #  X          vector of K forecasts

  # Output:
  #             Value of the BMA cdf evaluated at x

  if(length(sigma) == 1)
  {
    sigma = rep(sigma, length(a))
  }


  sum (w * pnorm (rep(Y,length(X)), a+b*X, sigma) )
}

