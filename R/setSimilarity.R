setSimilarity = function(candidateModules){
  
moduleCounts <- length(candidateModules)
HPI <- matrix(0,moduleCounts,moduleCounts)


for (i in 1:(moduleCounts-1)){
 module_i <- candidateModules[[i]]
 for (j in (i+1):moduleCounts){
  module_j <- candidateModules[[j]]
  if ((!(is.null(module_i)) && (!(is.null(module_j))))){
    HPI[i, j] <- length(intersect(module_i,module_j))/min(length(module_i),length(module_j))
  }
 }
}
#HPI[is.na(HPI)] <- 0
return(HPI)
}
