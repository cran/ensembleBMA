`dateCheck` <-
function (YYYYMMDDHH)
{
 if (!exists("chron")) library("chron")
 origin. <- c(month = 1, day = 1, year = 2000)
 YYYYMMDDHH <- sapply(YYYYMMDDHH, as.character)
 chk <- rep(TRUE, length(YYYYMMDDHH))
 l <- sapply(YYYYMMDDHH, nchar)
 if (any(L <- (l > 10 | l == 9 | l < 8))) chk[L] <- FALSE
 u <- unique(l[!L])
 dropHour <- if (length(u) > 1) FALSE else u == 8
 year <- as.numeric(sapply( YYYYMMDDHH[!L], substring, first = 1, last = 4))
 month <- as.numeric(sapply( YYYYMMDDHH[!L], substring, first = 5, last = 6))
 day <- as.numeric(sapply( YYYYMMDDHH[!L], substring, first = 7, last = 8))
 julianDate <- julian( month, day, year, origin. = origin.)
 ymdh <- julTOymdh( julianDate, origin = origin., dropHour = dropHour)
 chk[!L] <- ymdh == YYYYMMDDHH[!L]
 chk
}
