"[.ensembleData" <-
function (x, i, j) 
{ 
    ncolx <- ncol(x)
    matchCall <- match.call()
    matchCall[[1]] <- as.name("[.data.frame")
    if (missing(i)) matchCall$i <- 1:nrow(x)
    if (missing(j)) matchCall$j <- 1:ncolx
    matchCall$drop <- FALSE
    nForcs <- sum(sapply(names(x), 
                        function(nam) charmatch("forc", nam, nomatch = 0)))
    if (!missing(j)) {
      v <- (1:ncolx)
      names(v) <- names(x)
      j <- v[j]
      names(j) <- NULL
      if (any(j > nForcs)) 
        stop("column index must be confined to the forecasts")
      matchCall$j <- c(j, (nForcs+1):ncolx)
    }
    x <- eval(matchCall, parent.frame())
    x
}

