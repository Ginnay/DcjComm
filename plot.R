####################################################
############### Visualization output ###############
####################################################

########## Heatmap plot of number ##########

#Unify the name of cell types
#colna<-colnames(expr_l_r)
#colna<-gsub('BasalCells','Basal cells',colna)
#colna<-gsub('BCells','B cells',colna)
#colna<-gsub('TCells','T cells',colna)
#colna<-gsub('EndothelialCells','Endothelial cells',colna)
#colna<-gsub('LuminalEpithelialCells','Luminal cells',colna)
#colnames(expr_l_r)<-colna
library(stringr)
library(dplyr)
library(ggplot2)
source('visualization.R')
rows <- rep(rownames(expr_l_r), each = ncol(expr_l_r))
cols <- rep(colnames(expr_l_r), time = nrow(expr_l_r))
value <- c(t(expr_l_r))
result <- data.frame(rows, cols, value)

newdata<-subset(result, result$value > 0.5)
cluster_L<-str_split(newdata$cols, "-", simplify = T)[,2]
cluster_R<-str_split(newdata$cols, "-", simplify = T)[,1]
data <- data.frame(cluster_R,cluster_L)
d<-scSeqComm_heatmap_cardinality(data = data, title = "Number of ligand-receptor of cell communications")
d

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

