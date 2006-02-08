"[.ensembleData" <-
function (x, i, j) 
{ 
    ncolx <- ncol(x)
    matchCall <- match.call()
    matchCall[[1]] <- as.name("[.data.frame")
    if (missing(i)) matchCall$i <- 1:nrow(x)
    nForcs <- ensembleSize(x)
    if (!missing(j)) {
      v <- (1:ncolx)
      names(v) <- names(x)
      j <- v[j]
      names(j) <- NULL
      if (any(j > nForcs)) 
        stop("column index must be confined to the forecasts")
      if (any(duplicated(j))) 
        stop("repeated forecasts not allowed")
      if (nForcs < ncolx) {
        matchCall$j <- c(j, (nForcs+1):ncolx)
      }
     else  {
       matchCall$j <- j
      }
       
      nForcs <- length(j)
    }
    else matchCall$j <- 1:ncolx
    if (!missing(i)) {
      v <- (1:nrow(x))
      names(v) <- names(x)
      i <- v[i]
      names(i) <- NULL
      if (any(duplicated(i))) 
        stop("repeated entries not allowed")
    }
    matchCall$drop <- FALSE
    x <- eval(matchCall, parent.frame())
    attr(x, "ensembleSize") <- nForcs
    x
}

