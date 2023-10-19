#' Claculate the association score of receptor and pathways
#' @param TFscore  The calculated tf score.
#' @param tfweights The tfweights obtained obove.
#' @param method The computation method of score, eg "sum", "mean".
#' @param N_cores The cores used, default "1".
#' @param backend set the operating mode, default "doParallel"
#' @export
#'
Receptor_pathway_score <- function(TFscore, tfweights, method)
{
  #transcription factors
  TFs <- names(TFscore)[names(TFscore) %in% tfweights$tf]

  #retrieve TF score
  x <-lapply(TFscore[names(TFscore) %in% TFs], function(x){
    return(x$tfscore)
  })
  names(x) <- TFs
  tfscores <- data.frame(tfactor=names(x),tf_score=unlist(x))

  #transcription factor, its weight and score
  tf_values <- left_join(tfweights,tfscores, by=c("tf" = "tfactor"))

  tf_values <- tf_values[!is.na(tf_values$tf_score),]

  #computation of score
  if (method=="sum")
  {
    rp_score_ppr <- sum(tf_values$tf_score * tf_values$tf_PPR)
    rp_score <- sum(tf_values$tf_score)
  }
  if (method=="mean")
  {
    rp_score_ppr <- weighted.mean(x = tf_values$tf_score, w = tf_values$tf_PPR)
    rp_score <- mean(tf_values$tf_score)
  }
  #computation of R-Pathway score

  list_genes <- paste(unique(unlist(lapply(TFscore[names(TFscore) %in% TFs], function(x){
    return(x$genes)
  }))), collapse = ",")
  return(c(rp_score,rp_score_ppr,list_genes))
}

