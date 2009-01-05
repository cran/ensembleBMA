`[.ensembleBMAgamma0` <-
function(x,d) {
  d <- as.character(d) 
  if (any(!dateCheck(d))) stop("improperly specified dates")
  m <- match( d, names(x$nIter), nomatch = 0)
  if (any(!m)) stop("dates not matched in model")
  x$training$table <- x$training$table[d]
  x$prob0coefs <- x$prob0oefs[,,d,drop = FALSE]
  x$biasCoefs <- x$biasCoefs[,,d,drop = FALSE]
  x$varCoefs <- x$varCoefs[,d,drop = FALSE]
  x$weights <- x$weights[,d,drop = FALSE]
  x$nIter <- x$nIter[d]
  x
}

