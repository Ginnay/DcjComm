#' Compute the tfactors score in all cluster
#' @param regulatoryDB  The input regulatoryDB.
#' @param genes_score The TGpvalues obtained above.
#' @param cell_cluster The input cell_group.
#' @param cutoff Select genes after Bonferroni correction according to cutoff, default "0.05.
#' @export
#'
Tfactors_score_all_clusters <- function(regulatoryDB, genes_score, cell_cluster , cutoff){

  #condition of cutoff parameter
  if (cutoff < 0 | cutoff > 1)
    stop("cutoff must be a numeric value between 0 and 1")

  # name of each cell cluster
  cluster_names <- names(cell_cluster)

  # list that will contain, for each cluster, the values of transcription factors scores
  tf_scores <- list()

  # each cell cluster
  for(cl in cluster_names){
    # compute the transcription factor score
    genes_pvalue <- genes_score[[cl]]$pvalue
    names(genes_pvalue) <- row.names(genes_score[[cl]])
    tf_scores[[cl]] <- Tfactors_score(regulatoryDB, genes_pvalue , cutoff)#
  }
  return(tf_scores)

}
