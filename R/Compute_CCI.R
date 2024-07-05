#' Classify genes into clusters
#' @param data  The input expression matrix.
#' @param cell_group The clusters obtained by cell cluster.
#' @export
#'
Compute_CCI = function(data = data, Label = Label, cell_group = cell_group, LR_pairs_DB = LR_pairs_DB,
                       TF_reg_DB = TF_TG, R_TF_association = TF_PPR,  LR_db = LR_db,
                       N_cores = 2, DEmethod = "wilcoxon", n_bootstrap = 4, backend = "doParallel",
                       other_inter_scores = NULL, cell_reference = NULL, Org = Org,
                       use.type="median",probs = 0.9,method="weighted"){



  ########## intracellular communication ##########

  library(doParallel)
  library(dplyr)
  library(methods)
  library(pheatmap)
  library(stringr)
  library(graphics)
  library(psych)
  library(clusterProfiler)
  library(stats)
  library(magrittr)
  library(utils)
  library(circlize)
  library(gridBase)
  library(grid)
  library(ComplexHeatmap)
  library(ggalluvial)
  library(ggplot2)
  library(RColorBrewer)
  library(jsonlite)
  library(networkD3)
  library(reshape2)
  library(tidyr)
  library(enrichplot)
  library(ggrepel)
  library(ggridges)
  library(DOSE)

  gene_expr<-as.data.frame(gene_expr)


  message("Check input data...")
  message("**** scRNA-seq gene expression matrix having ")
  message(" - ", nrow(gene_expr), " genes")
  message(" - ", ncol(gene_expr), " cells")
  message(" - ", length(cell_group), " cell groups/clusters",
          " (  ", paste0(names(cell_group), sep = "  "), ")")
  genes <- rownames(gene_expr)
  message("**** Trascriptional regulatory networks database having ",
          length(TF_reg_DB), " transcription factors and ", length(unique(unlist(TF_reg_DB))),
          " genes")
  regulatoryDB <- lapply(X = TF_reg_DB, FUN = function(X,
                                                       gene_IDs) {
    tgenes <- X[X %in% genes]
    return(tgenes)
  }, gene_IDs = genes)
  regulatoryDB <- regulatoryDB[!(lapply(regulatoryDB, length) ==
                                   0)]
  Tf <- unique(names(regulatoryDB))
  tg <- unique(unlist(regulatoryDB))
  message("**** Considering only transcription factors regulating genes present in the input scRNA-seq data, the trascriptional regulatory networks is reduced to ",
          length(Tf), " transcription factors and ", length(tg),
          " genes")
  message("\nAnalyzing intracellular communications...")
  message("**** Compute TF activity")
  message(" - Identify realiable target genes")
  tictoc::tic()
  classification_genes <- assign_genes_all_clusters(exprMatr = gene_expr,
                                                    cell_clusters = cell_group, H0_cells = cell_reference,
                                                    N_cores = N_cores, backend = backend)
  tictoc::toc()
  gc()
  message(" - Identify differentially expressed target genes")
  tictoc::tic()
  source('Pvalue_target_genes_all_cluster.R')
  TGpvalue <- Pvalue_target_genes_all_cluster(exprMatr = gene_expr,
                                              cell_cluster = cell_group, target_genes = genes, H0_cells = cell_reference,
                                              unreliable_genes_all_clusters = classification_genes,
                                              N_cores = N_cores, backend = backend, DEmethod = DEmethod)
  tictoc::toc()
  rm(classification_genes)
  gc()
  message(" - Compute TF activity score")
  tictoc::tic()
  TF_activity <- Tfactors_score_all_clusters(regulatoryDB = regulatoryDB,
                                             genes_score = TGpvalue, cell_cluster = cell_group, cutoff = 0.05)
  tictoc::toc()
  message("**** Load TF PPR scores")
  PPR_scores <- R_TF_association
  message("**** Compute intracellular signaling evidence (i.e. score S_intra)")
  tictoc::tic()
  S_intra_all_clusters <- S_intra_all_clusters(R_pathway_TF_weight = PPR_scores,
                                               TFscore = TF_activity, cell_cluster = cell_group, method = "mean")
  tictoc::toc()
  S_intra_all_clusters <- cbind(S_intra_all_clusters, up_genes = NA,
                                down_genes = NA)
  for (i in 1:nrow(S_intra_all_clusters)) {
    list_genes <- strsplit(S_intra_all_clusters[i, "list_genes"],
                           ",")[[1]]
    x <- TGpvalue[[S_intra_all_clusters[i, "cluster"]]][list_genes,
    ]
    up_genes <- row.names(x)[which(x$sign == "up")]
    down_genes <- row.names(x)[which(x$sign == "down")]
    S_intra_all_clusters[i, "up_genes"] <- paste(up_genes,
                                                 collapse = ",")
    S_intra_all_clusters[i, "down_genes"] <- paste(down_genes,
                                                   collapse = ",")
  }
  S <- subset(S_intra_all_clusters, S_intra_all_clusters$receptor %in%
                LR_db$receptor)
  comm_results <- left_join(S, LR_db, by = "receptor")
  comm_results <- comm_results[, c("row", "cluster", "S_intra")]
  S_intra <- matrix(data = comm_results[, 3], nrow = length(unique(comm_results[,
                                                                                1])), ncol = length(unique(comm_results[, 2])), byrow = TRUE,
                    dimnames = list(unique(comm_results[, 1]), unique(comm_results[,
                                                                                   2])))
  S_intra <- as.data.frame(S_intra)
  row <- gsub("_", "-", rownames(S_intra))
  rownames(S_intra) <- row
  S_intra <- S_intra[rowSums(S_intra[]) > 0, ]

  ########## intercellular communication ##########
  B <- c(1:dim(gene_expr)[2])
  C <- paste(B[1:dim(gene_expr)[2]],Label[1:dim(gene_expr)[2]],sep = "_")
  colnames(gene_expr) <- C
  gene_expr <- as.data.frame(gene_expr)

  library(stringr)
  col<-colnames(gene_expr)
  col1<-str_split(C, "_", simplify = T)[,2]
  unique(col1)

  mt <- CreateNichConObject(data=gene_expr, min.feature = 3,
                            names.field = 2,
                            names.delim = "_",
                            source = "TPM",
                            scale.factor = 10^6,
                            Org = Org,
                            project = "Microenvironment")


  object = mt

  complex_tmp <- LR_db$receptor[grep(",",LR_db$receptor)] %>% unique()
  tmp_complex_symbol <- LR_db$receptor[grep(",",LR_db$receptor)] %>% unique() %>% str_split(",") %>% unlist %>% unique()
  all.gene.needed <- unique(as.character(c(LR_db$ligand, LR_db$receptor, R_TF_association$tf, TF_reg_DB$TF, TF_reg_DB$TG,tmp_complex_symbol)))

  my_Expr <- object@data$withoutlog
  colnames(my_Expr) <- as.character(object@meta.data$celltype)
  my_Expr[1:4,1:4]
  detect_gene <- rownames(my_Expr)

  expr_set <- my_Expr[intersect(detect_gene, all.gene.needed),]
  detect_gene <- rownames(expr_set)
  cell_type = unique(colnames(expr_set))
  expr.fc <- object@data$withoutlog[detect_gene,]
  colnames(expr.fc) <- colnames(expr_set)

  rm(list=c("object"))

  complex_matrix <- matrix(ncol = length(colnames(expr_set)))
  complex_matrix <- as.data.frame(complex_matrix)
  colnames(complex_matrix) <- colnames(expr_set)
  myrownames <- c()

  complex <- complex_tmp
  if(length(complex)>0){
    for(i in 1:length(complex)){
      i_tmp = strsplit(complex[i], ',')
      # print(i_tmp)
      if( sum(i_tmp[[1]] %in% detect_gene) == length(i_tmp[[1]]) ){
        tmp_df <- expr_set[i_tmp[[1]],]
        tmp_mean <- colMeans(tmp_df)
        tmp_index <- unique(unlist(apply(tmp_df, 1,function(x) {which(x==0)})))
        tmp_mean[tmp_index] <- 0
        complex_matrix <- rbind(complex_matrix, tmp_mean)
        myrownames <- c(myrownames, complex[i])
      }
    }

    complex_matrix <- complex_matrix[-1,]

    if(nrow(complex_matrix) > 0){
      rownames(complex_matrix) <- myrownames
      expr_set <- rbind(expr_set, complex_matrix)
    }
  }

  expr_set <- expr_set[apply(expr_set, 1, function(x){sum(x!=0)})>0,]
  detect_gene <- rownames(expr_set)
  # expr_set[1:4,1:4]

  expr_mean <- matrix(nrow = nrow(expr_set), ncol = length(cell_type))
  myColnames <- c()
  for (i in 1:length(cell_type)) {
    myCell <- cell_type[i]
    myMatrix <- expr_set[,colnames(expr_set)==myCell,drop=F]
    if(use.type=="mean"){
      myMatrix_mean <- as.numeric(apply(myMatrix, 1, mean))
    }else if(use.type=="median"){
      quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
        quantile(x, probs = probs,names=FALSE)
      }))
      mean.tmp <- rowMeans(myMatrix)
      mean.tmp[which(quantil.tmp==0)]<-0
      myMatrix_mean <- mean.tmp
    }
    expr_mean[,i] <- myMatrix_mean
    myColnames <- c(myColnames, myCell)
    # print(myCell)
  }
  expr_mean <- data.frame(expr_mean)
  colnames(expr_mean) <- myColnames
  rownames(expr_mean) <- rownames(expr_set)

  expr_mean <- expr_mean[apply(expr_mean, 1, function(x){sum(x!=0)})>0,]
  detect_gene <- rownames(expr_mean)

  if(use.type=="median"){
    fc.list <- mylog2foldChange.diy(inData = expr.fc, cell.type = cell_type, method="median", probs = probs)
  }else{
    fc.list <- mylog2foldChange.diy(inData = expr.fc, cell.type = cell_type, method="mean", probs = probs)
  }

  l_r_inter <- unique(LR_db[,2:3])
  colnames(l_r_inter)<-c("Ligand_Symbol","Receptor_Symbol")

  # softmax for ligand
  ligand_symbol <- unique(LR_db$ligand)
  softmax_ligand <- expr_mean[intersect(ligand_symbol, detect_gene),]
  colnames(softmax_ligand) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_ligand)
  softmax_ligand <- do.call(rbind,lapply(1:nrow(softmax_ligand), function(i){
    softmax_ligand[i,]/rowCounts[i]
  }))

  # softmax for receptor
  receptor_symbol <- unique(LR_db$receptor)
  softmax_receptor <- expr_mean[intersect(receptor_symbol, detect_gene),]
  colnames(softmax_receptor) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_receptor)
  softmax_receptor <- do.call(rbind,lapply(1:nrow(softmax_receptor), function(i){
    softmax_receptor[i,]/rowCounts[i]
  }))

  #  l-r in cell type level
  l_r_inter <- unique(LR_db[,2:3])
  colnames(l_r_inter)<-c("Ligand_Symbol","Receptor_Symbol")
  expr_l_r <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)^2) ##A->A,A->B,A->C,,,,C->C
  expr_l_r <- as.data.frame(expr_l_r)
  rownames(expr_l_r) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  myColnames <- character()
  for (i in cell_type) {
    for (j in cell_type) {
      myColnames <- c(myColnames, paste(i,j,sep = "-"))
    }
  }
  colnames(expr_l_r) <- myColnames

  ######
  expr_sender <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type))
  rownames(expr_sender) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames( expr_sender) <- cell_type

  expr_receiver <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type))
  rownames(expr_receiver) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(expr_receiver) <- cell_type

  expr_tf <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)) ##A->A,A->B,A->C,,,,C->C
  expr_tf  <- as.data.frame(expr_tf)
  rownames(expr_tf ) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(expr_tf)<-cell_type

  sender_val_weighted <- 0
  receiver_val_weighted <- 0

  for (n in 1:nrow(l_r_inter)) {
    sender_tmp <- l_r_inter[n,1]
    receiver_tmp <- l_r_inter[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    #print(n)
    for (i in cell_type) {
      for (j in cell_type) {
        myColnames <- c(myColnames, paste(i,j,sep = "-"))
        val_tmp = 0
        #val_tmp3 = 0

        if( sum(l_r_inter[n,] %in% detect_gene)==2 ){
          sender_val <- expr_mean[sender_tmp,i]
          receiver_val <- expr_mean[receiver_tmp,j]

          if(sender_val>0 & receiver_val >0){
            sender_val_weighted <- softmax_ligand[sender_tmp, i]
            receiver_val_weighted <- softmax_receptor[receiver_tmp, j]
            val_tmp <- 100*(sender_val_weighted^2 + receiver_val_weighted^2) * S_intra
            #val_tmp <- sqrt((sender_val_weighted + receiver_val_weighted) * S_intra)
          }else{
            val_tmp = 0
          }
        }else{
          val_tmp = 0
        }
        #####
        col_index_tmp <- paste(i,j,sep = "-")
        expr_sender[n,i] <- sender_val_weighted
        expr_receiver[n,j] <- receiver_val_weighted
        expr_l_r[n,col_index_tmp] <- val_tmp
      }
    }
  }

  expr_l_r_log2 <- log2(expr_l_r+1)
  expr_l_r_log2_scale <- (expr_l_r_log2-min(expr_l_r_log2))/(max(expr_l_r_log2)-min(expr_l_r_log2))

  expr_l_r <- list(expr_l_r_log2_scale = expr_l_r_log2_scale)
  return(expr_l_r)
}
