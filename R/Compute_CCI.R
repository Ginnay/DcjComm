#' Compute the cell-cell interaction
#' @param gene_expr  The input expression matrix.
#' @param TF_PPR The signal transductions from receptor to transcription factors.
#' @param TF_TG The transcription factor-target regulatory interactions.
#' @param S_intra The calculated results of intracellular score.
#' @param row_index_final The ligand-receptor association.
#' @param Org The input species information.
#' @export
#'
Compute_CCI = function(gene_expr = gene_expr,TF_PPR=TF_PPR,TF_TG=TF_TG, S_intra=S_intra,row_index_final=row_index_final,
                       Org = Org, use.type="median",probs = 0.9,method="weighted"){

  mt <- CreateNichConObject(data=gene_expr, min.feature = 3,
                          names.field = 2,
                          names.delim = "_",
                          source = "TPM",
                          scale.factor = 10^6,
                          Org = Org,
                          project = "Microenvironment")


  object = mt
  library(stringr)
  complex_tmp <- row_index_final$receptor[grep(",",row_index_final$receptor)] %>% unique()
  tmp_complex_symbol <- row_index_final$receptor[grep(",",row_index_final$receptor)] %>% unique() %>% str_split(",") %>% unlist %>% unique()
  all.gene.needed <- unique(as.character(c(row_index_final$ligand, row_index_final$receptor, TF_PPR$tf, TF_TG$TF, TF_TG$TG,tmp_complex_symbol)))

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

  l_r_inter <- unique(row_index_final[,2:3])
  colnames(l_r_inter)<-c("Ligand_Symbol","Receptor_Symbol")

  # softmax for ligand
  ligand_symbol <- unique(row_index_final$ligand)
  softmax_ligand <- expr_mean[intersect(ligand_symbol, detect_gene),]
  colnames(softmax_ligand) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_ligand)
  softmax_ligand <- do.call(rbind,lapply(1:nrow(softmax_ligand), function(i){
    softmax_ligand[i,]/rowCounts[i]
  }))

  # softmax for receptor
  receptor_symbol <- unique(row_index_final$receptor)
  softmax_receptor <- expr_mean[intersect(receptor_symbol, detect_gene),]
  colnames(softmax_receptor) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_receptor)
  softmax_receptor <- do.call(rbind,lapply(1:nrow(softmax_receptor), function(i){
    softmax_receptor[i,]/rowCounts[i]
  }))

  #  l-r in cell type level
  l_r_inter <- unique(row_index_final[,2:3])
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
  expr_Sinter<-expr_l_r
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
        val_Sinter = 0
        #val_tmp3 = 0

        if( sum(l_r_inter[n,] %in% detect_gene)==2 ){
          sender_val <- expr_mean[sender_tmp,i]
          receiver_val <- expr_mean[receiver_tmp,j]

          if(sender_val>0 & receiver_val >0){
            sender_val_weighted <- softmax_ligand[sender_tmp, i]
            receiver_val_weighted <- softmax_receptor[receiver_tmp, j]
            val_tmp <- 100*(sender_val_weighted^2 + receiver_val_weighted^2) * S_intra
            val_Sinter <- 100*(sender_val_weighted^2 + receiver_val_weighted^2)
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
        expr_Sinter[n,col_index_tmp] <- val_Sinter
      }
    }
  }

  expr_l_r_log2 <- log2(expr_l_r+1)
  expr_l_r_log2_scale <- (expr_l_r_log2-min(expr_l_r_log2))/(max(expr_l_r_log2)-min(expr_l_r_log2))

  expr_Sinter_log2 <- log2(expr_Sinter+1)
  expr_Sinter_log2_scale <- (expr_Sinter_log2-min(expr_Sinter_log2))/(max(expr_Sinter_log2)-min(expr_Sinter_log2))

  expr_l_r <- list(expr_l_r_log2_scale = expr_l_r_log2_scale,expr_Sinter=expr_Sinter_log2_scale)
  return(expr_l_r)
}
