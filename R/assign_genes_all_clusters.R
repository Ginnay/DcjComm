#' Classify genes into clusters
#' @param exprMatr  The input expression matrix.
#' @param cell_clusters The clusters obtained by cell cluster.
#' @param H0_cells The specific list of reference cells, default "NULL".
#' @param N_cores The cores used, default "1".
#' @param backend set the operating mode, default "doParallel"
#' @export
#'
assign_genes_all_clusters <- function(exprMatr, cell_clusters, H0_cells = NULL, N_cores = 1, backend = "doParallel")
{
  # if no specific list of reference cells is specified, the whole matrix is taken as reference
  if(is.null(H0_cells)){
    H0_cells = colnames(exprMatr)
  }

  #object where results will be stored
  classification_genes <- list()

  #names of clusters
  cluster_names <- names(cell_clusters)

  if (backend == "doParallel") doParallel::registerDoParallel(cores = N_cores)
  if (backend == "doMC") doMC::registerDoMC(cores = N_cores)

  #message(paste("... using",foreach::getDoParWorkers(),"cores with",foreach::getDoParName(),"..."))

  classification_genes <- foreach(cl_id=1:length(cluster_names), .packages = "Matrix") %dopar% {

    cl <- cluster_names[[cl_id]]
    reliable_genes_cluster <- rownames(exprMatr) [(Matrix::rowSums(exprMatr[,cell_clusters[[cl]]]>1) / length(cell_clusters[[cl]])) > 0.25]

    #gene classified as reliable in the reference
    if(is.list(H0_cells)){
      H0_cells_cl <- H0_cells[[cl]]
    }else{
      H0_cells_cl <- setdiff(H0_cells,cell_clusters[[cl]])
    }


    reliable_genes_ref <- rownames(exprMatr) [(Matrix::rowSums(exprMatr[,H0_cells_cl]>1) / length(H0_cells_cl)) > 0.25]

    reliable_genes <- union(reliable_genes_cluster, reliable_genes_ref)

    #unreliable genes
    unreliable_genes <- setdiff(rownames(exprMatr),reliable_genes)

    #results
    list(reliable = reliable_genes, unreliable= unreliable_genes)

  }

  gc(verbose = FALSE)
  names(classification_genes) <- cluster_names

  if (backend == "doParallel") doParallel::stopImplicitCluster()

  return(classification_genes)
}


