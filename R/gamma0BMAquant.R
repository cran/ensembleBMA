"gamma0BMAquant" <-
function(alpha, WEIGHTS, PROB0, MEAN, VAR)
{
 # if the probability of zero is greater than the desired alpha
 # then the quantile is zero

  if (sum(WEIGHTS*PROB0) > alpha) return(0)

# Initialize: Find lower and upper bounds
 
  lower <- 0
  upper <- max(MEAN+6*sqrt(VAR))

  if (gamma0BMAcdf(lower, WEIGHTS, PROB0, MEAN, VAR, 0) > alpha) return(NA)
  if (gamma0BMAcdf(upper, WEIGHTS, PROB0, MEAN, VAR, 0) < alpha) return(NA)

  z <- uniroot(gamma0BMAcdf, lower = lower, upper = upper,
               WEIGHTS=WEIGHTS, PROB0=PROB0, MEAN=MEAN, VAR=VAR, offset = alpha)

# print(c(alpha, z$root,abs(z$f.root)))

  z$root
}

