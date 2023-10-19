#' Compute the ligand_receptor communications in all cluster
#' @param exprMatr  The input expression matrix.
#' @param cell_clusters The clusters obtained by cell cluster.
#' @param ligands The input target genes.
#' @param receptors The target genes score.
#' @export
#'
Ligand_receptor_all_cluster_score <- function(exprMatr, cell_cluster, ligands, receptors, H0_genes = NULL, verbose = FALSE, debug_result = FALSE, bootstrap = FALSE, bootstrap_size = 1){

  # name of each cell cluster
  cluster_names <- names(cell_cluster)

  # list that will contain, for each cluster, the values of ligands and receptors scores
  ligands_receptors_S_score <- list()

  # identify row indices of ligands and receptors
  lig_id <- which(row.names(exprMatr) %in% ligands)
  rec_id <- which(row.names(exprMatr) %in% receptors)

  # if no specific list of "background" genes is specified, H0 will be built using all the genes
  if(is.null(H0_genes)){
    #H0_genes = rownames(exprMatr)
    H0_genes <- 1:nrow(exprMatr)
  }

  # if array of gene names, the row indices are identified
  if(!is.numeric(H0_genes))
  {
    H0_genes <- which(row.names(exprMatr) %in% H0_genes)
  }


  # each cell cluster
  for(cl in cluster_names){

    cell_id <- which(colnames(exprMatr) %in% cell_cluster[[cl]])

    if (bootstrap)
    {
      # bootstrapping sampling of columns
      cell_id <- sample(cell_id, size = round(length(cell_id)*bootstrap_size), replace = TRUE)
      print("bootstrapping")
    }

    # compute the ligands and receptors S score for the subset of the count matrix related to the current cell cluster
    ligands_receptors_S_score[[cl]] <- compute_score_S_ligand_receptor1 (exprMatr = exprMatr,
                                                                        cells = cell_id,
                                                                        ligands = lig_id,
                                                                        receptors = rec_id,
                                                                        H0_genes = H0_genes,
                                                                        verbose = verbose, debug_result = debug_result)
  }

  return(ligands_receptors_S_score)
}
