#' Select highly variable genes
#' @param params The input params.
#' @export
#'
Selectgene = function(params){

  library(RSoptSC)
  gene_expression_threshold <- params$threshold
  n_features <- params$number
  data <- params$data
  filtered_data <- SelectData(data, gene_expression_threshold, n_features)
  X <- filtered_data$M_variable
  X <- log2(X+1)

  Topgene <- list(X = X)
  return(Topgene)
}
