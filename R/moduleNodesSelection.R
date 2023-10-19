#' Assign genes to different clusters
#' @param moduleparams  The input params.
#' @export
#'
moduleNodesSelection = function( moduleparams){

 U_final <- moduleparams$U_final
 xita <- moduleparams$xita

 N <- dim(U_final)[1]
 K <- dim(U_final)[2]

 candidateModules <- vector('list',K)
 moduleSignal <- numeric(K)

 H_mean <- apply(U_final,2,mean)
 H_std <- apply(U_final,2,sd)

 for (k in 1:K){

  candidateModules[[k]] <- which(U_final[, k] > H_mean[k] + xita*H_std[k])
  moduleSignal[k] <- mean(U_final[candidateModules[[k]], k])
 }

 HPI <- setSimilarity( candidateModules )
 modulesFinal <- candidateModules
 HPI[is.na(HPI)] <- 0
 mHPI <- sum(t(sum(HPI)))/sum(sum(HPI!=0))


 module_result <- list(modulesFinal=modulesFinal, mHPI=mHPI)

 return(module_result)

}
