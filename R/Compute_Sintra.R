#' Compute the intercellular interaction
#' @param gene_expr  The input expression matrix.
#' @param cell_group The input cell_group.
#' @param LR_pairs_DB The input ligand-receptor interactions.
#' @param TF_reg_DB The transcription factor-target regulatory interactions.
#' @param R_TF_association The input PPR_scores.
#' @param LR_db The input ligand-receptor interactions.
#' @export
#'
Compute_Sintra = function(gene_expr = gene_expr, cell_group = cell_group,LR_pairs_DB = LR_pairs_DB,
                          TF_reg_DB = TF_TG, R_TF_association = TF_PPR, LR_db=row_index_final,
                           N_cores = 2,DEmethod = "wilcoxon",n_bootstrap = 4,backend = "doParallel",
                           other_inter_scores = NULL,cell_reference = NULL){
  library(doParallel)
  library(dplyr)
  library(methods)
  library(pheatmap)
  library(stringr)
  message("Check input data...")

  message("**** scRNA-seq gene expression matrix having ")
  message(" - ", nrow(gene_expr), " genes")
  message(" - ", ncol(gene_expr), " cells")
  message(" - ", length(cell_group), " cell groups/clusters", " (  ",paste0(names(cell_group), sep = "  " ), ")")


  # get genes IDs
  genes <- rownames(gene_expr)

  message("**** Trascriptional regulatory networks database having ", length(TF_reg_DB), " transcription factors and ", length(unique(unlist(TF_reg_DB))), " genes")

  # remove the entries in the transcriptional regulatory networks database corresponding to genes not present in the given scRNA-seq data
  regulatoryDB <- lapply(X=TF_reg_DB,
                         FUN=function(X, gene_IDs){
                           tgenes <- X[X %in% genes]
                           return(tgenes)
                         },
                         gene_IDs = genes)
  # remove transcription factors for which all the regulated genes are not present in the given scRNA-seq data
  regulatoryDB <- regulatoryDB[!(lapply(regulatoryDB,length)==0)]

  # get remaining transcription factors (Tf) and target genes (tg)
  Tf <- unique(names(regulatoryDB))
  tg <- unique(unlist(regulatoryDB))

  message("**** Considering only transcription factors regulating genes present in the input scRNA-seq data, the trascriptional regulatory networks is reduced to ", length(Tf), " transcription factors and ", length(tg), " genes")

  message("\nAnalyzing intracellular communications...")

  message("**** Compute TF activity")

  # Identification of realible genes to be used to measure the effect of intracellular signaling
  message(" - Identify realiable target genes")

  tictoc::tic()
  classification_genes <-  assign_genes_all_clusters(exprMatr = gene_expr,
                                                             cell_clusters = cell_group,
                                                             H0_cells = cell_reference,
                                                             N_cores = N_cores,
                                                             backend = backend)
  tictoc::toc()
  gc()

  # Computation of DE genes for each cluster. Genes in each cluster are compared with a referencce group of cell
  message(" - Identify differentially expressed target genes")
  tictoc::tic()
  TGpvalue <- Pvalue_target_genes_all_cluster(exprMatr = gene_expr,
                                                      cell_cluster = cell_group,
                                                      target_genes = genes,
                                                      H0_cells = cell_reference,
                                                      unreliable_genes_all_clusters = classification_genes,
                                                      N_cores = N_cores,
                                                      backend = backend,
                                                      DEmethod = DEmethod)
  tictoc::toc()
  rm(classification_genes); gc()

  # Computation of TFscore
  message(" - Compute TF activity score")
  tictoc::tic()
  TF_activity  <- Tfactors_score_all_clusters(regulatoryDB = regulatoryDB,
                                                      genes_score = TGpvalue,
                                                      cell_cluster = cell_group,
                                                      cutoff = 0.05)
  tictoc::toc()

  #computation of R-Pathway score in each cluster
  message("**** Load TF PPR scores")

  PPR_scores <- R_TF_association

  message("**** Compute intracellular signaling evidence (i.e. score S_intra)")
  tictoc::tic()
  S_intra_all_clusters <- S_intra_all_clusters(R_pathway_TF_weight = PPR_scores,
                                                             TFscore = TF_activity,
                                                             cell_cluster = cell_group,
                                                             method="mean")
  tictoc::toc()

  #add info about sign of DE analysis previously computed
  S_intra_all_clusters <- cbind(S_intra_all_clusters, up_genes = NA, down_genes = NA)

  for (i in 1:nrow(S_intra_all_clusters))
  {
    # info about DE genes and cluster
    list_genes <- strsplit(S_intra_all_clusters[i, "list_genes"],",")[[1]]
    x <- TGpvalue[[S_intra_all_clusters[i, "cluster"]]][list_genes,]

    # save list of up and down-regulated genes
    up_genes <- row.names(x)[which(x$sign == "up")]
    down_genes <- row.names(x)[which(x$sign == "down")]
    S_intra_all_clusters[i,"up_genes"] <- paste(up_genes,collapse = ",")
    S_intra_all_clusters[i,"down_genes"] <- paste(down_genes,collapse = ",")
    #if (length(list_genes) != length(c(up_genes, down_genes))) print("help")
  }

  S <- subset(S_intra_all_clusters,S_intra_all_clusters$receptor%in%LR_db$receptor)

  comm_results <- left_join(S,LR_db, by="receptor")
  comm_results<-comm_results[,c("row","cluster","S_intra")]

  S_intra <- matrix(data = comm_results[,3], nrow = length(unique(comm_results[,1])), ncol = length(unique(comm_results[,2])), byrow = TRUE, dimnames = list(unique(comm_results[,1]), unique(comm_results[,2])))
  S_intra <-as.data.frame(S_intra)
  row<-gsub('_','-',rownames(S_intra))
  rownames(S_intra)<-row

  S_intra <- S_intra[rowSums(S_intra[])>0,]

  Intra <- list(S_intra = S_intra, TF_Score = TF_activity)
  return(Intra)
}
