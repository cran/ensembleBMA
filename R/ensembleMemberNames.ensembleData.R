"ensembleMemberNames.ensembleData" <-
function (x) 
{ 
 i <- charmatch("obs", dimnames(x)[[2]]) - 1
 (dimnames(x)[[2]])[1:i]
}

