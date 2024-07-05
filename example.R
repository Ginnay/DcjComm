####################################################
#### PART 0: Data collection and preprocessing #####
####################################################

# Load DcjComm package and example data
setwd("E:/DcjComm/DcjComm/R") ##Set working directory
library(DcjComm) ##load the DcjComm package
library('Matrix')
library('scRMD')
library(doParallel)
library(dplyr)

data("lr_mouse") ##ligand-receptor pairs databases (L-R)
LR_pairs_DB = lr_mouse[,c("ligand","receptor")]
data("TF_PPRmouse") ##Receptor-Transcription factor a-priori association(TF_PPR)
data("TF_TGmousedb") ##Transcriptional regulatory networks (TF_TG)
data("lst_scrna") ##Load the input data
data <- lst_scrna$mat
data <- as.matrix(data)
data("cell_group") ##Load the cell group
data("true_labs") ##The true cell clster
data("W") ##Load the symmetric weight matrix calculated by graph regularization

#optional--elimate batch effect
#library(sva)
#data <- ComBat(dat=as.matirx(data), batch=batch)

# Select highly variable genes
threshold <- 0.06 ##parameter of filtering the low-quality cells
number <- 2000  ##highly variable genes
params <- list(data = data, threshold = threshold, number = number) ##input parameters
Selectgene(params) -> Topgene ##select highly variable genes
X <- Topgene[["X"]]


####################################################
################ PART 1: Cell cluster ##############
####################################################

params <- list(X = X, W = W, k1 = 100, true_labs = true_labs) ##Set the input parameters of cell cluster
NMFbased(params) -> OUT_result ##Cell cluster using the NMF-based joint learning method
Label<-OUT_result[["Label"]] ##The cell labels afer cell clustering

#Annotate cell clusters
#Label[Label == 1] = "BCells" ##'Cd79a', 'Cd79b', 'Cd74', 'Cd19'# B --cluster1
#Label[Label == 2] = "Macrophages" ##"Cd74","Cd14","Csf1r"# Macro --cluster2
#Label[Label == 3] = "LuminalEpithelialCells" ##"Krt8","Krt18","Krt19"# Luminal --cluster3
#Label[Label == 4] = "EndothelialCells" ##'Pecam1', 'Nrp1', 'Kdr','Oit3' #Endo --cluster4
#Label[Label == 6] = "BasalCells"  ##"Krt14","Krt5","Krt17")) # Basal --cluster6
#Label[Label == 5] = "Fibroblasts" ##number of cell types
#Label[Label == 7] = "TCells" ##number of cell types


####################################################
########## PART 2: Cell-cell communication #########
####################################################

# Select marker genes by Seurat
#library(tidyverse)
#library(Matrix)
#library(Seurat)
#data<-log2(data+1)
#pbmc <- CreateSeuratObject(counts = data, project = "scrna", min.cells = 0, min.features = 0)
#group<-Label
#pbmc@meta.data[["RNA_snn_res.0.06"]]<-group
#pbmc@meta.data[["seurat_clusters"]]<-group
#Idents(object = pbmc) <- "RNA_snn_res.0.06"
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#row<-unique(pbmc.markers$gene,pbmc.markers$gene %in% rownames(data))
#gene_expr <- data[row,] ##Infer cell-cell communications by marker genes

LR_db <- lr_mouse
TF_reg_DB <- TF_TGmousedb
R_TF_association <- TF_PPRmouse
Org <- "Mus musculus" #choose the species source of gene, eg "Homo sapiens", "Mus musculus"

# Infer cell-cell communications
CCI <- Compute_CCI(
  data = gene_expr, Label = Label, cell_group = cell_group, LR_pairs_DB = LR_pairs_DB,
  TF_reg_DB = TF_reg_DB, R_TF_association = R_TF_association, LR_db = LR_db,
  N_cores = 2, DEmethod = "wilcoxon", n_bootstrap = 4, backend = "doParallel",
  other_inter_scores = NULL, cell_reference = NULL, Org = Org,
  use.type="median",probs = 0.9,method="weighted"
)
expr_l_r <-CCI$expr_l_r_log2_scale ##The results of cell-cell communications

###################################################
######### PART 3: Determination of modules #########
####################################################

# Select module
U_final<-OUT_result[["U_final"]] ##The matrix to detect modulesS
xita = 2 ##The input parameter
moduleparams <- list(U_final=U_final, xita=xita) ##The input information
moduleNodesSelection(moduleparams) -> module_result
S_final<-OUT_result[["S_final"]] ##The matrix to calculate the importance of modules


