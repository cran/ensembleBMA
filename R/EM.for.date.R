"EM.for.date" = function(date, date.list, X, Y, eps = 1e-005, maxiter=1000, start.w=NULL, start.sigma=NULL, const.var = TRUE, num.training.days = 25, lead=2, reg.adjust = TRUE, min.CRPS = TRUE)
#EM.for.date fits a model to forecast a given date based on a 25-day training sample
#with a 1-day gap (so the forecast is 48 hours - if we give it day index of 27, it trains on days 1 through 25)
#X is the matrix of forecasts, Y is the vector of observations, date is the date you want to forecast for,
#and date.list is the list of the dates corresponding to each observation
#the specific format for the dates in date and date.list does not matter, so long as they are numeric, and
#a difference of 1 indicates a difference of 1 day
#this function discards from the training set any observations for which we are missing one or more of the member forecasts
{
  foo=(1:length(Y))
  baz=rep(0, length(Y))
  for(i in 1:(dim(X)[2]))
  {
    baz=baz+is.na(X[,i])
  }
  bar=foo[((date-(lead + num.training.days))<date.list) & (date.list<(date-lead+1)) & (baz==0)]
  output=EM.normals(X[bar,1:(dim(X)[2])], Y[bar], eps, maxiter,
                    start.w, start.sigma, const.var, reg.adjust, min.CRPS)
  output
}
