`[.ensembleBMAnormal` <-
function (x, d) 
{
    d <- as.character(d)
    if (any(!dateCheck(d))) 
        stop("improperly specified dates")
    m <- match(d, names(x$nIter), nomatch = 0)
    if (any(!m)) 
        stop("dates not matched in model")
    x$training$table <- x$training$table[d]
    x$biasCoefs <- x$biasCoefs[, , d, drop = FALSE]
    x$sd <- if (is.null(dim(x$sd))) x$sd[d] else x$sd[,d]
    x$weights <- x$weights[, d, drop = FALSE]
    x$nIter <- x$nIter[d]
    attr(x, "call") <- list(attr(x, "call"), match.call())
    x
}
