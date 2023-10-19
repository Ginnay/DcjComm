#' Compute the pvalue of target genes in all cluster
#' @param exprMatr  The input expression matrix.
#' @param cell_clusters The clusters obtained by cell cluster.
#' @param target_genes The input target genes.
#' @export
#'
Pvalue_target_genes_all_cluster <- function(exprMatr, cell_cluster, target_genes, DEmethod = "MAST", H0_cells = NULL, unreliable_genes_all_clusters, N_cores = 1, backend = "doParallel")
{

  # if no specific list of reference cells is specified, the wilcoxon will be computed using all the cells
  if(is.null(H0_cells)){
    H0_cells = colnames(exprMatr)
  }

  # name of each cell cluster
  cluster_names <- names(cell_cluster)

  # list that will contain, for each cluster, the values of target genes scores
  tg_scores <- list()

  gc(verbose = FALSE)

  if (DEmethod != "MAST")
  {
    if (backend == "doParallel") doParallel::registerDoParallel(cores = N_cores)
    if (backend == "doMC") doMC::registerDoMC(cores = N_cores)
    message(paste("... using",foreach::getDoParWorkers(),"cores with",foreach::getDoParName(),"..."))
  }else{
    message(paste("... using parallelism of MAST:",N_cores, "cores ..."))
    if (backend == "doParallel") doParallel::registerDoParallel(cores = 1)
    if (backend == "doMC") doMC::registerDoMC(cores = 1)
  }
  tg_scores <- foreach(cl_id=1:length(cluster_names), .packages = "Matrix") %dopar% {

    cl <- cluster_names[[cl_id]]

    if(is.list(H0_cells)){
      H0_cells_cl <- H0_cells[[cl]]
    }else{
      H0_cells_cl <- setdiff(H0_cells,cell_cluster[[cl]])
    }

    # perform DE analysis
    source('Pvalue_target_genes.R')
    Pvalue_target_genes(exprMatr = exprMatr[target_genes, cell_cluster[[cl]]],
                                ref_exprMatr = exprMatr[target_genes, H0_cells_cl],
                                target_genes = target_genes,
                                unreliable_genes = unreliable_genes_all_clusters[[cl]]$unreliable,
                                DEmethod = DEmethod,
                                N_cores = N_cores)
  }
  names(tg_scores) <- cluster_names

  if (backend == "doParallel") doParallel::stopImplicitCluster()

  return(tg_scores)
}
