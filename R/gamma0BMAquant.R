"gamma0BMAquant" <-
function(alpha, WEIGHTS, POP, MEAN, VAR)
{
 # if the probability of zero is greater than the desired alpha
 # then the quantile is zero

  if (sum(WEIGHTS*POP) > alpha) return(0)

# Initialize: Find lower and upper bounds
 
  lower <- 0
  upper <- max(MEAN+6*sqrt(VAR))

  if (gamma0BMAcdf(lower, 0, WEIGHTS, POP, MEAN, VAR) > alpha) return(NA)
  if (gamma0BMAcdf(upper, 0, WEIGHTS, POP, MEAN, VAR) < alpha) return(NA)

  z <- uniroot(gamma0BMAcdf, lower = lower, upper = upper,
               offset=alpha, WEIGHTS=WEIGHTS, POP=POP, MEAN=MEAN, VAR=VAR)

# print(c(alpha, z$root,abs(z$f.root)))

  z$root
}

