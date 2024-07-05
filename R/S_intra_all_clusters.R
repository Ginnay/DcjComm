#' Compute the cell-cell interaction
#' @param R_pathway_TF_weight  The input PPR_scores.
#' @param TFscore The obtained TF_activity.
#' @param cell_cluster The above constructed cell_group.
#' @export
#'
S_intra_all_clusters <- function(R_pathway_TF_weight, TFscore, cell_cluster, method=c("sum","mean"))
{

  Rpathway_score <- data.frame(cluster = character(0), receptor = character(0) , pathway = character(0),
                               rpathway_score = numeric(0), S_intra = numeric(0), list_genes = character(0))

  #all pairs R-Pathway
  rp_pair <- R_pathway_TF_weight[!duplicated(R_pathway_TF_weight[,c("receptor","pathway")]),c("receptor","pathway")]

  # name of each cell cluster
  cluster_names <- names(cell_cluster)

  #receptor-pathway score depends on clusters
  for (cl in cluster_names)
  {
    #for each R-Pathway pair
    for (j in 1:nrow(rp_pair))
    {
      #trascription factors and its weights associated to the current R-Pathway pair
      tfweights <- R_pathway_TF_weight[R_pathway_TF_weight$receptor == rp_pair[j,]$receptor & R_pathway_TF_weight$pathway == rp_pair[j,]$pathway,c("tf","tf_PPR")]
      if (!any(names(TFscore[[cl]]) %in% tfweights$tf)) next

      #computation of score
      scores <- Receptor_pathway_score(TFscore[[cl]],tfweights,method=method)

      #save results
      Rpathway_score <- rbind(Rpathway_score , data.frame(cluster = cl, receptor = rp_pair[j,"receptor"], pathway = rp_pair[j,"pathway"],
                                                          rpathway_score = as.numeric(scores[1]), S_intra = as.numeric(scores[2]),
                                                          list_genes = as.character(scores[3])))
      Rpathway_score$list_genes <- as.character(Rpathway_score$list_genes)
    }
  }

  metadata_info <- setdiff(colnames(R_pathway_TF_weight), c("receptor","tf_PPR","tf","pathway"))
  if(length(metadata_info)!=0)
  {
    metadata <- R_pathway_TF_weight[,!(colnames(R_pathway_TF_weight) %in% c("receptor","tf_PPR","tf"))]
    metadata <- metadata[!duplicated(metadata),]
    Rpathway_score <- left_join(Rpathway_score, metadata, by= "pathway")
  }


  return(Rpathway_score)
}
