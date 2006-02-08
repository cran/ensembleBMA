"cdfBMA.fitBMAgamma0" <-
function(fit, ensembleData, values, dates = NULL, popData = NULL, ...) 
{

 M <- matchEnsembleMembers(ensembleData, fit)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

 inverseLogit <- function(x) {
# logit function safeguared against underflow and overflow
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  1/(1+exp(-x))
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  if (x >= log(.Machine$double.eps)) {
                    x <- exp(x)
                    x/(1+x)
                  }
                  else exp(x)
                }
                else 0
             }
            }

# inverseLogit <- function(x) exp(x)/(1 + exp(x))

 if (!is.null(popData) && !is.null(dim(popData))) {
   if (length(dim(popData)) == 2) {
     popData <- list(popData)
   }
   else {
     popData <- apply(popData, 3, list)
   }
 }

 nObs <- ensembleNobs(ensembleData)

 ensDates <- ensembleDates(ensembleData)
 if (!is.null(dates) && length(dates) > 1 && is.null(ensDates)) 
   stop("date ambiguity")
 
 L <- 1:nObs

 if (!is.null(dates) && !is.null(ensDates)) {
   M <- as.logical(match(as.character(ensDates), dates, nomatch=0))
   if (!any(M)) stop("dates not matched in data")
   L <- L[M]
 }

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData), as.character(values)) 

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 if (!any(is.na(WEIGHTS <- fit$weights)))  {

    for (i in L) {
    
       f <- ensembleData[i,]
     
       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f
        
       fTrans <- sapply(f, fit$transformation)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       if (is.null(popData)) {
         PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs,
                              2,sum), inverseLogit)
       }
       else {
         popi <- rbind(lapply( popData, function(x,i) x[i,], i = i))
         PROB0 <- sapply(apply(rbind( 1, fTrans, popi) * fit$prob0coefs,
                                    2,sum), inverseLogit)
       }

       CDF[i,] <- sapply( sapply( values, fit$transformation), 
                          gamma0BMAcdf, 
              WEIGHTS = WEIGHTS, POP = POP, MEAN = MEAN, VAR = VAR) 
    }

 }

 CDF[ L, , drop=FALSE]
}

