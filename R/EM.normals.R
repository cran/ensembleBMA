
"EM.normals"= function(X, Y, eps = 1e-005, maxiter=1000, start.w=NULL, start.sigma=NULL, const.var = TRUE, reg.adjust = TRUE, min.CRPS = TRUE)
# function for EM algorithm for mixture of ensemble-member-centered normals
# Modification of Adrian's modification of Fadoua's EM code for SLP - Adrian 9/11/03, 10/21/03
# McLean 12/01/03 - modifying to use for 8-member wind speed (should be generalizeable to any mixture of normals)
# McLean 05/30/04 - incorporating Veronica's changes to have a single variance parameter for all models
# McLean 06/18/04 - fixed matrix transpose errors in calculating z

# Inputs:
#  X           matrix of ensemble members. This is an n by K matrix, where there are
#              n observations to be used in the fitting, and K ensemble members
#  Y           n-vector of observations
#  eps         stopping criterion
#  maxiter     maximum number of EM iterations allowed
#  start.w     initial values for the weights (optional)
#  start.sigma initial value for the sd

{
  n=length(Y)
  k=dim(X)[2]
  lik = 0
  z=matrix(ncol=k, nrow=n)
  w.new=rep(0, k)
  v.new=rep(0, k)

  start.v=(start.sigma)^2

  # set initial value for v, either as total variance or given starting value
  if(is.null(start.sigma))
  {
     v=var(Y)
  }
  else
  {
    v=start.v
  }

  # set intial weights, either as equal weights or given starting values
  if(is.null(start.w))
  {
    w=rep(1/k, k)
  }
  else
  {
    w=start.w
  }

  error=rep(1, 3)
  niter=0

  if(reg.adjust == TRUE)
  {
    #bias-correction terms (linear regression)
    B=rep(0, k)
    A=rep(0, k)
    for(i in 1:k)
    {
      B[i]=var(X[,i], Y, na.rm=TRUE)/var(X[,i], na.rm=TRUE)
      A[i]=mean(Y)-B[i]*mean(X[,i])
    }
  }
  else {
    B=rep(1, k)
    A=rep(0, k)
  }

  if(const.var == TRUE) {
  #main EM algorithm
  while((max(abs(error)) > eps) && (niter < maxiter))
  {
    sumz=0
    z.old=z
    #weighted sum of distribution functions for each member (sumz is a vector, one value for each observation)
    for(i in 1:k)
    {
      sumz=sumz+w[i]*dnorm(Y, mean=A[i]+B[i]*X[,i], sd=sqrt(v))
    }
    #matrix of probabilities for each observation coming from each member (latent variables)
    z=t(w*t(dnorm(Y, mean=t(A+B*t(X)), sd=sqrt(v))))/sumz    

    #calculate new weights and variance based on latent variables
    for(i in 1:k)
    {
      w.new[i]=sum(z[,i])/n
    }
    v.new=sum(z*((Y-t(A+B*t(X)))^2))/sum(z)
    
    # Compute log-likelihood (corrected by Adrian on 10/21/03)
    lik.old=lik
    lik=sum(log(sumz))

    #calculate change from last iteration to this one
    error[1]=max(abs(w.new-w))
    error[2]=max(abs(log(v.new/v)))
    if(niter==0)
    {
      error[3]=1
    }
    else
    {
      error[3]=max(abs(z-z.old))
    }
    v=v.new
    w=w.new
    niter <- niter + 1
    
  }
  }
 else {

  v = rep(var(Y),k)

  while((max(abs(error)) > eps) && (niter < maxiter))
  {
    sumz=0
    #weighted sum of distribution functions for each member (sumz is a vector, one value for each observation)

    #matrix of probabilities for each observation coming from each member (latent variables)
      z.old=z
	for(i in 1:k)
	{
		value <- dnorm(Y, mean=A[i]+B[i]*X[,i], sd=sqrt(v[i]))* w[i]
		sumz <- sumz + value
		z[,i] <- value
	}

		for( i in 1:k)
		{
			z[,i] <- z[,i]/sumz
	}


    #calculate new weights and variance based on latent variables
    for(i in 1:k)
    {
      w.new[i]=sum(z[,i])/n
	v.new[i]=sum(z[,i]*(Y-A[i]-B[i]*X[,i])^2)/sum(z[,i])
    }
    
    # Compute log-likelihood (corrected by Adrian on 10/21/03)
    lik.old=lik
    lik=sum(log(sumz))

    #calculate change from last iteration to this one
    error[1]=max(abs(w.new-w))
    error[2]=max(abs(log(v.new/v)))
    if(niter==0)
    {
      error[3]=1
      error[4]=1
    }
    else
    {
      error[3]=max(abs(z-z.old))
      error[4]=max(abs(lik-lik.old))
    }

    v=v.new
    w=w.new
    niter <- niter + 1

  }
  }

  if(min.CRPS==TRUE)
  {
    if(const.var==TRUE)
    {
      sigma.CRPS=sqrt(v[1])
      CRPS.optim=function(sigma)
      {
        CRPS(A,B,rep(sigma,k),w,X,Y)
      }
      v=rep((optimize(CRPS.optim,interval=c(0, 2*sigma.CRPS))$minimum)^2,k)
    }
    else
    {
      sigma.CRPS=sqrt(v)
      CRPS.optim=function(sigma)
      {
        CRPS(A,B,sigma,w,X,Y)
      }
      v=(optim(sigma.CRPS, CRPS.optim)$par)^2
    }
  }

  list(loglik=lik, a=A, b=B, w=w, sigma=sqrt(v), z=z, niter=niter)
}
