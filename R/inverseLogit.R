"inverseLogit" <-
function(x) {

## inverse logit function exp(x) / (1 + exp(x))
## safeguared against underflow and overflow

              if (is.na(x)) return(NA)
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  x <- exp(-x)
                  1/(1+x)
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  x <- exp(x)
                  x/(1+x)
                }
                else 0
             }
            }

