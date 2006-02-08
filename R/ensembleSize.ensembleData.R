"ensembleSize.ensembleData" <-
function(x) 
{ 
 charmatch("obs", dimnames(x)[[2]]) - 1
}

