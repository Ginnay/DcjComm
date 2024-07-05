#' Compute the tfactors score
#' @param regulatoryDB  The input regulatoryDB.
#' @param genes_score The TGpvalues obtained above.
#' @param cutoff Select genes after Bonferroni correction according to cutoff, default "0.05.
#' @export
#'
Tfactors_score <- function(regulatoryDB, genes_score , cutoff)
{
  #list of all genes
  genes <- names(genes_score)

  #all genes with score < cutoff after Bonferroni correction
  Tup <- names(genes_score[genes_score < (cutoff/length(genes))])

  #compute transcription factor score as 1-pvalue of Fisher test
  tfscores <- lapply(regulatoryDB, function(x) { Ltf <- intersect(x,genes)

  #perform Fisher test
  tf_score <- fisher_score(Tup,Ltf,genes)

  #DE genes
  DE_tgenes <- intersect(Tup,Ltf)

  #results
  res <- list(tfscore = tf_score, genes = DE_tgenes)

  return(res)})

  return(tfscores)
}
