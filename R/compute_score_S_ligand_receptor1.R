#' Compute the pvalue of target genes in all cluster
#' @param exprMatr  The input expression matrix.
#' @param cells The input cell_id.
#' @param ligands The input lig_id.
#' @param cells The input cell_id.
#' @param receptors The input rec_id.
#' @param  H0_genes The input array of gene names.
#' @export
#'
compute_score_S_ligand_receptor1 <- function(exprMatr, cells, ligands, receptors, H0_genes, verbose = TRUE, debug_result = FALSE) {

  # if verbose, print some infos
  if(verbose == TRUE){
    message("Num. cells: ", ncol(exprMatr))
    message("Num. ligands: ", length(ligands))
    message("Num. receptors: ", length(receptors))
    message("Num. H0 genes: ", length(H0_genes))

    if(length(H0_genes) == nrow(exprMatr)){
      message("H0 will include all the gene in the matrix")
    }
    else{
      message("H0 will include a subset of the genes in the matrix")
    }
  }

  # compute ligands and receptors expression levels (as mean of their normalized counts across cells) and compute mean and sd of gaussian approximation
  library(scSeqComm)
 ### res_cpp <- BigBootstrap(pBigMat = exprMatr@address,
 ###                        bootstrap_ind = cells,
 ###                         H0_genes = H0_genes,
 ###                        ligands = ligands,
 ###                         receptors = receptors)

  # # compute ligands and receptors expression levels (as mean of their normalized counts across cells)
   lig_avg_expr <- Matrix::rowMeans(exprMatr[ligands, ])
   rec_avg_expr <- Matrix::rowMeans(exprMatr[receptors, ])
  #
  # # compute mean and sd of gaussian approximation
   #H0_mean <- Matrix::mean(exprMatr[H0_genes,])
   #H0_sd <- sd(exprMatr[H0_genes,]) / sqrt(N_COL)

  score_L <- pnorm(lig_avg_expr, mean = 0, sd = 1)
  score_R <- pnorm(rec_avg_expr, mean = 0, sd = 1)

  # imposing zero score when average expression levels is zero
  if (length(which(lig_avg_expr ==0 ))>0) score_L[which(lig_avg_expr == 0)] <- 0
  if (length(which(rec_avg_expr ==0 ))>0) score_R[which(rec_avg_expr == 0)] <- 0

  names(score_L) <- row.names(exprMatr)[ligands]
  names(score_R) <- row.names(exprMatr)[receptors]

  results <- list(ligand_score = score_L, receptor_score = score_R)

  # add mean and sd of H0 to function return value
  if(debug_result == TRUE){
    results$H0 <- list(H0_mean = H0_mean, H0_sd = H0_sd)
  }

  return(results)

}
