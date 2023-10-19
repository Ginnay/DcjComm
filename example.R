####################################################
########## PART 0: Load packages and data ##########
####################################################

########## Load DcjComm package ##########
setwd("/data/dingqian/DcjComm/DcjComm/R")
library(DcjComm)

# load example data and database
data("lr_mouse") #ligand-receptor pairs databases (L-R)
data("TF_PPRmouse") # Receptor-Transcription factor a-priori association(TF_PPR)
data("TF_TGmousedb") # Transcriptional regulatory networks (TF_TG)

data("lst_scrna") #input data
data("cell_group")
data("true_labs")
data("W")
data <- lst_scrna$mat
data <- as.matrix(data)

####################################################
########## PART 1: Select highly variable genes ##########
####################################################
library('Matrix')
library('scRMD')
threshold <- 0.06 #Parameter of filtering the low-quality cells
number <- 2000  #highly variable genes

params <- list(data = data, threshold = threshold, number = number)
Selectgene(params) -> Topgene
X <- Topgene[["X"]]

####################################################
########## PART 2: Dimension Reduction and cell cluster ##########
####################################################
k1 <- 100
k2 <- length(unique(true_labs))

params <- list(X = X, W = W, k1 = k1, k2 = k2)
DrjccLG(params) -> OUT_result
F_final<-OUT_result[["F_final"]]

Label=0
for (e in 1:dim(F_final)[2]){
  v=F_final[,e]
  ma=max(v)
  s=which(v==ma)
  Label[e]=s
}
Label<-as.matrix(Label)

ARI<-CalculateARI(true_labs, Label)

####################################################
########## PART 3: Detection cell-cell communications ##########
####################################################
########## Select DEGs for detecting cell-cell communication ##########
library(tidyverse)
library(Matrix)
library(Seurat)

data<-log2(data+1)
pbmc <- CreateSeuratObject(counts = data, project = "scrna", min.cells = 0, min.features = 0)

########## Assign cell cluster according to marker genes and number of cell types ##########
Label[Label == 1] = "BCells" ##'Cd79a', 'Cd79b', 'Cd74', 'Cd19'# B --cluster1
Label[Label == 2] = "Macrophages" ##"Cd74","Cd14","Csf1r"# Macro --cluster2
Label[Label == 3] = "LuminalEpithelialCells" ##"Krt8","Krt18","Krt19"# Luminal --cluster3
Label[Label == 4] = "EndothelialCells" ##'Pecam1', 'Nrp1', 'Kdr','Oit3' #Endo --cluster4
Label[Label == 6] = "BasalCells"  ##"Krt14","Krt5","Krt17")) # Basal --cluster6
Label[Label == 5] = "Fibroblasts" ##number of cell types
Label[Label == 7] = "TCells" ##number of cell types

group<-Label
pbmc@meta.data[["RNA_snn_res.0.06"]]<-group
pbmc@meta.data[["seurat_clusters"]]<-group
Idents(object = pbmc) <- "RNA_snn_res.0.06"

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
row<-unique(pbmc.markers$gene,pbmc.markers$gene %in% rownames(data))
gene_expr <- data[row,]

########## intracellular communication ##########
LR_pairs_DB <- lr_mouse[,c("ligand","receptor")]
TF_TG <- TF_TGmousedb
TF_PPR <- TF_PPRmouse
R_TF_association <- TF_PPR


library(doParallel)
library(dplyr)
library(methods)
library(pheatmap)
library(stringr)
gene_expr<-as.data.frame(gene_expr)
Intra_score <- Compute_Sintra(
  gene_expr = gene_expr, cell_group = cell_group,LR_pairs_DB = LR_pairs_DB,
  TF_reg_DB = TF_TG, R_TF_association = TF_PPR, LR_db=lr_mouse
)

S_intra<-Intra_score[["S_intra"]]
TF_activity<-Intra_score[["TF_Score"]]
gc()

########## intercellular communication ##########
B <- c(1:4481)
C <- paste(B[1:4481],Label[1:4481],sep = "_")
colnames(gene_expr) <- C
gene_expr <- as.data.frame(gene_expr)

library(stringr)
col<-colnames(gene_expr)
col1<-str_split(C, "_", simplify = T)[,2]
unique(col1)

library(graphics)
library(stringr)
library(methods)
library(dplyr)
library(psych)
library(clusterProfiler)
library(stats)
library(magrittr)
library(utils)
library(circlize)
library(gridBase)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
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
Org = "Mus musculus"
CCI <- Compute_CCI(
  gene_expr = gene_expr,TF_PPR=TF_PPR,TF_TG=TF_TG, S_intra=S_intra,row_index_final=lr_mouse,
  Org = Org
)
expr_l_r <-CCI$expr_l_r_log2_scale
expr_Sinter<-CCI[["expr_Sinter"]]
####################################################
########## PART 4: Determination of modules ##########
####################################################
###module
U_final<-OUT_result[["U_final"]]
xita = 2
moduleparams <- list(U_final=U_final, xita=xita)
moduleNodesSelection(moduleparams) -> module_result

S_final<-OUT_result[["S_final"]]

####################################################
########## PART 5: Visualization output ##########
####################################################

########## Heatmap plot of number ##########
colna<-colnames(expr_l_r)
colna<-gsub('BasalCells','Basal cells',colna)
colna<-gsub('BCells','B cells',colna)
colna<-gsub('TCells','T cells',colna)
colna<-gsub('EndothelialCells','Endothelial cells',colna)
colna<-gsub('LuminalEpithelialCells','Luminal cells',colna)

colnames(expr_l_r)<-colna
colnames(expr_Sinter)<-colna

rows <- rep(rownames(expr_l_r), each = ncol(expr_l_r))
cols <- rep(colnames(expr_l_r), time = nrow(expr_l_r))
value <- c(t(expr_l_r))
result <- data.frame(rows, cols, value)
library(stringr)
newdata<-subset(result, result$value > 0.5)
cluster_L<-str_split(newdata$cols, "-", simplify = T)[,2]
cluster_R<-str_split(newdata$cols, "-", simplify = T)[,1]
data <- data.frame(cluster_R,cluster_L)

#path<-"F:/DrjCCI/code"
#setwd(path)
source('visualization.R')

library(dplyr)
library(ggplot2)
d<-scSeqComm_heatmap_cardinality(data = data, title = "Number of ligand-receptor of cell communications")
d


########## Heatmap plot of interaction ##########
y = "cluster_L"
x = "cluster_R"
limit_fill = NULL
cell_types = NULL
#select variables from data
data <- data[,c(x,y)]
#count the number of observations of x-y pair
df <- data %>% count(data[[x]],data[[y]])

if (is.null(cell_types)) cell_types = unique(c(data[[x]],data[[y]]))
res <- gtools::permutations(n=length(cell_types),r=2,v=cell_types,repeats.allowed=T)
colnames(df) <- c(x,y,"cardinality")
colnames(res) <- c(x,y)
df <- left_join(as.data.frame(res[,1:2]), as.data.frame(df))
df[is.na(df)]=0
count <- matrix(data = df[,3], nrow = 7, ncol = 7, byrow = TRUE, dimnames = list(unique(df[,1]), unique(df[,2])))
groupSize <- as.numeric(table(Label))

color.use<-NULL
library(ComplexHeatmap)
library(CellChat)
mat<-count

if (is.null(color.use)) {
  color.use <- scPalette(ncol(mat))
}
names(color.use) <- colnames(mat)

color.heatmap<-"Reds"
color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))
row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

if (sum(abs(mat) > 0) == 1) {
  color.heatmap.use = c("white", color.heatmap.use)
} else {
  mat[mat == 0] <- NA
}

legend.name<-"Number of interactions"
font.size.title<-10
title.name<-"Number of interactions"
font.size<-8
cluster.rows<-FALSE

ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
              bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = cluster.rows,cluster_columns = cluster.rows,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
              # width = unit(width, "cm"), height = unit(height, "cm"),
              column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA, #at = colorbar.break,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)

ht1

########### Sankey plot of Ligand-receptor-tf interaction ##########
result<-result[which(result[,3]>0.5),]

library(stringr)
ligand<-str_split(result$rows, "-", simplify = T)[,1]
receptor<-str_split(result$rows, "-", simplify = T)[,2]
result<-cbind(ligand,receptor,result$value)
result<-as.data.frame(result)

library(dplyr)
comm<- left_join(result,TF_PPR, by="receptor")

comm<- na.omit(comm)
tmp<-comm[,c("ligand","receptor","tf","tf_PPR")]
colnames(tmp)<-c("Ligand","Receptor","tf","value")
tmp<-unique(tmp)

#important ligands Tgfa, Ngf, Col5a2, Il11, Col4a2, Jag1, Col18a1  Hspg2)
fData1 = tmp[grepl('Tgfa',tmp$Ligand),]
fData2 = tmp[grepl('Ngf',tmp$Ligand),]
fData3 = tmp[grepl('Col5a2',tmp$Ligand),]
fData4 = tmp[grepl('Tl11',tmp$Ligand),]
fData5 = tmp[grepl('Col4a2',tmp$Ligand),]
fData6 = tmp[grepl('Jag1',tmp$Ligand),]
fData7 = tmp[grepl('Col18a1',tmp$Ligand),]
fData8 = tmp[grepl('Hspg2',tmp$Ligand),]

#important receptor Gpc1, Procr, Fzd7, Itga5, Ldlr, Tlr2, Lrp6, Ephb1, and Tfrc
fData9 = tmp[grepl('Gpc1',tmp$Receptor),]
fData10 = tmp[grepl('Procr',tmp$Receptor),]
fData11 = tmp[grepl('Fzd7',tmp$Ligand),]
fData12 = tmp[grepl('Itga5',tmp$Receptor),]
fData13 = tmp[grepl('Ldlr',tmp$Ligand),]
fData14 = tmp[grepl('Tlr2',tmp$Receptor),]
fData15 = tmp[grepl('Lrp6',tmp$Ligand),]
fData16 = tmp[grepl('Ephb1',tmp$Receptor),]
fData17 = tmp[grepl('Tfrc',tmp$Receptor),]

fData = rbind(fData1,fData2,fData3,fData4,fData5,fData6,fData7,
              fData8,fData9,fData10,fData11,fData12,fData13,
              fData14,fData15,fData16,fData17)
fData = unique(fData)

tmp<-fData

tmp1<-tmp[order(tmp$value,decreasing = T),][1:50,]
tmp.df<-tmp1

mycol.vector = c('#5d62b5','#29c3be','#f2726f','#62b58f','#bc95df', '#67cdf2', '#ffc533', '#5d62b5', '#29c3be')
elments.num <-  tmp.df %>% unlist %>% unique %>% length()
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

df<-tmp1
library(ggplot2)
tmp.df
axes=1:3
mycol = mycol.vector.list[1:elments.num]
font.size = 2
boder.col="white"
set_alpha = 0.8
subdf <- df
colnames(subdf) <- c("Ligand_symbol", "Receptor_symbol", "TF", "value")
subdf <- as.data.frame(subdf)
options(stringsAsFactors = FALSE)
diy_stratum <- ggalluvial::StatStratum

p <- ggplot(as.data.frame(subdf),
            aes(y = value,
                axis1 = Ligand_symbol,
                axis2 = Receptor_symbol,
                axis3 = TF)) +
  scale_fill_manual(values = mycol) +
  ggalluvial::geom_flow(stat = "alluvium",width = 1/8,aes(fill = Ligand_symbol), alpha=set_alpha) +
  ggalluvial::geom_stratum(width = 1/8, reverse = T,alpha = .9, size =0.001) +
  geom_text(stat = diy_stratum,
            #family = "A",
            size = 2.6, face = "bold", color = "black", aes(label = after_stat(stratum)),
            reverse = T) +
  scale_x_continuous(breaks = 1:3, labels = c("Ligand", "Receptor", "TF")) +

  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(
          #family = "A",
          size = 4, face = "bold", color = "black")) +

  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(text = element_text(
    #family = "A",
    size = 10, face = "bold", color = "black"))+
  ggtitle("")
p

########### Circos plot ##########
library(gridBase)
library(circlize)
cell_color <- data.frame(color=c('#e31a1c','#1f78b4','#e31a1c','#1f78b4',
                                 'peachpuff','lightgrey','lightpink'
), stringsAsFactors = FALSE)
rownames(cell_color) <- c("B cells","Macrophages","T cells","Basal cells",
                          "Luminal cells","Endothelial cells","Fibroblasts")

ViewInterCircos(object = expr_l_r, font = 2,
                cellColor = cell_color,
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 2,
                order.vector=c("B cells","Macrophages","T cells","Basal cells",
                               "Luminal cells","Endothelial cells","Fibroblasts"),
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)

########### ridge plot ##########
library(ggplot2)
library(ggridges)
library(viridis)
g<-ggplot(fdata, aes(x = tfscore , y = fenzu , fill = fenzu)) +
  geom_density_ridges_gradient(scale=3)+
  theme_ridges()+
  xlab(NULL) + ylab(NULL) +  DOSE::theme_dose()+
  theme(text = element_text(
    size = 10, face = "bold", color = "black",family="Times")) + theme_bw()+
  theme(axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.text.x = element_text(size = 9, face = "bold", color = "black"))
g




