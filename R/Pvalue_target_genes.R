#' Compute the cell-cell interaction
#' @param gene_expr  The input expression matrix.
#' @param TF_PPR The signal transductions from receptor to transcription factors.
#' @param TF_TG The transcription factor-target regulatory interactions.
#' @param S_intra The calculated results of intracellular score.
#' @param row_index_final The ligand-receptor association.
#' @param Org The input species information.
#' @export
#'
Pvalue_target_genes <- function(exprMatr, ref_exprMatr, target_genes, unreliable_genes, DEmethod = "MAST", N_cores = 1)
{

  results <- data.frame(pvalue = rep(1,times = length(target_genes)), sign = NA)
  # #array(1, dim = length(target_genes))
  rownames(results) <- target_genes

  test_genes <- setdiff(target_genes, unreliable_genes)

  if (DEmethod == "MAST")
  {
    matr <- cbind(as.matrix(exprMatr[test_genes,]), as.matrix(ref_exprMatr[test_genes,]))
    colnames(matr) <- 1:ncol(matr)
    grp <- c(rep(1,ncol(exprMatr)),rep(2,ncol(ref_exprMatr)))
    cdr <- scale(colMeans(matr > 0))

    # construct object
    sca <- suppressMessages(MAST::FromMatrix(exprsArray = matr,
                                             cData = data.frame(wellKey = colnames(matr),
                                                                grp = grp, cdr = cdr)))
    #parallelism options
    options(mc.cores = N_cores)

    #fit the model
    zlmdata <- suppressMessages(MAST::zlm(~cdr + grp, sca, parallel = T))
    mast <- suppressMessages(MAST::lrTest(zlmdata, "grp"))

    #sign of DE
    for (gg in test_genes){
      results[gg,"pvalue"] <- mast[gg, "hurdle", "Pr(>Chisq)"]
      results[gg,"sign"] <- ifelse(mean(exprMatr[gg,]) > mean(ref_exprMatr[gg,]), "up", "down")
    }
  }


  if (DEmethod == "wilcoxon")
  {

    CHUNCK_SIZE <- 1000
    test_genes_list <- split(test_genes, ceiling(seq_along(test_genes) / CHUNCK_SIZE))
    for (ch in test_genes_list){
      if(length(ch)==1)
      {
        tmp_1 <- t(as.matrix(exprMatr[ch,]))
        row.names(tmp_1) <- ch
        tmp_2 <- t(as.matrix(ref_exprMatr[ch,]))
        row.names(tmp_2) <- ch
      }else{
        tmp_1 <- as.matrix(exprMatr[ch,])
        tmp_2 <- as.matrix(ref_exprMatr[ch,])
      }

      # wilcoxon test and sign of DE
      for (gg in ch){
        results[gg,"pvalue"] <- wilcox.test(x=tmp_1[gg,],y=tmp_2[gg,],exact=FALSE)$p.value
        results[gg,"sign"] <- ifelse((mean(tmp_1[gg,])>mean(tmp_2[gg,])), "up", "down")
      }
      rm(tmp_1, tmp_2); gc(verbose = FALSE)
    }

  }

  gc(verbose = FALSE)

  return(results)
}

