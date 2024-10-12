
#---------Expression data from esophageal squamous cell carcinoma patients---------------
### Esophageal Cancer vs Paratumor (CGC Congress 2023- www.cgc2023.com) ###
#GSE161533  #Samples (84 total) 
#28 Paratumor and 28 Tumor specimens
#GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

#Load Packages:
library(Biobase)
library(limma)
library(GEOquery)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

#Load previous data to speed up:
setwd("G:\\OOF\\Bioinformatics\\Projects\\2nd Congress\\Microarray\\Data")
lapply(X = c("DEGs.RData","EXDATA.RData","HEATMAPOO.RData"),load,.GlobalEnv)

#--------------------- LOAD DATA -------------------------
setwd("G:\\OOF\\Bioinformatics\\Projects\\2nd Congress\\Microarray\\Data")

series <- "GSE39144" 
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "dataq/")
platform<- "GPL570"
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#Export Expression Matrix (Expression data) :
ex <- exprs(gset)
dim(ex)


#--------------------- SAMPLE SELECTION ---------------------
#Sample Selection:
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111111111111111111",
               "1111110000000000000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# X = Normal  1 = paratumor 0 = tumor

# filter out excluded samples (marked as "X")
sel <- which(sml != "X") #here you can choose which sample you want to exclude
sml <- sml[sel]
gset <- gset[ ,sel]
#Rewrite expression for selected samples:
##in this stage we define ex without norma samples
ex <- exprs(gset)
dim(ex)

#Log2 transformation needed?
max(ex)
min(ex)
#Yes or No???
## if YES:
ex<- log2(ex+1) #The data is scaled!
# Apply ITTTTT!!!
exprs(gset) <- ex

#for plots:
gr <- factor(c(rep("Paratumor",28), rep("Tumor", 28)))

#Save Data:
setwd("E:/Projects/Mammad Crush CGC/RData/")
save(data, ex, gr, file = "EXDATA.RData")

#--------------------- QC PLOTS ---------------------
setwd("E:\\Projects\\CGC2023/CGC Mammad_Esophageal Cancer/plots/")

### Draw BoxPlot:
pdf("boxplot.pdf")
boxplot(ex)
dev.off()


### Draw PCA:
#Method 1
exm <- ex-rowMeans(ex) #rowMeans is a vector 
pc <- prcomp(exm)
pcr <-data.frame(pc$r[,1:3], gr)


png("PCA_rowMeans_N vs C.png",width = 2000, height = 1500,res = 300,units = "px")

ggplot(pcr,aes(PC1,PC2,color=gr)) + ggtitle("PCA rowMeans")+
  geom_label(size=2.2, fontface = "bold",label= gr)+
  geom_vline(xintercept=0, col="black",linetype = "longdash") + 
  geom_hline(yintercept=0, col="black",linetype = "longdash") + 
  stat_ellipse(aes(fill=gr),type = "norm",geom = "polygon",
               alpha = 0.09,level = 0.6)+
  theme_bw()

dev.off()

#Method 2
ex.scale <- t(scale(t((ex)), scale=T,  center = F))
pc <- prcomp(ex.scale)
pcr <- data.frame(pc$r[,1:3], Group=gr)


png("PCA_scale_N vs C.png",width = 2000, height = 1500,res = 300,units = "px")

ggplot(pcr,aes(PC1,PC2,color=gr)) + ggtitle("PCA scaled")+
  geom_label(size=2.2, fontface = "bold",label= gr)+
  geom_vline(xintercept=0, col="black",linetype = "longdash") + 
  geom_hline(yintercept=0, col="black",linetype = "longdash") + 
  stat_ellipse(aes(fill=gr),type = "norm",geom = "polygon",
               alpha = 0.09,level = 0.6)+
  theme_bw()


dev.off()


#--------------------- Find DEGs ---------------------
setwd("E:\\Projects\\Mammad Crush CGC\\results")
getwd()
## Make Design Matrix:
# assign samples to groups and set up design matrix
gs <- factor(sml) #   0 = Cancer , 1 = Normal
groups <- make.names(c("Tumor","Paratumor")) #the order of the names should be by Levels(gs) order
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

#Limma starts:
fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-") 
cts #check: "Tumor-Paratumor" / "Cancer-Normal" / "Experiment-Control"

cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust.method ="fdr", sort.by="logFC", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","logFC","adj.P.Val","P.Value","t","B","ID"))

head(tT)

#remove non gene probes:
tT <- tT[-which(tT$Gene.symbol == ""), ]

#remove /// in gene symbols and select first name:
tT$Gene.symbol <- sub("///.*","",tT$Gene.symbol)
#remove /// in gene IDs and select first ID:
tT$Gene.ID <- sub("///.*","",tT$Gene.ID)

#write.table(tT, file="total results.txt", quote=F, sep="\t")

# Total DEGs(up and down): for main heatmap

#DEGs:#  IDs
degs <- subset(tT, abs(logFC) > 2 & adj.P.Val < 0.05)

###  multiple probe sets mapping to the same gene: A need for aggregation (Thanks to Dr. Bahari) ###

nums <- table(degs[,"Gene.symbol"])
degs.ag <- data.frame(matrix(0,nrow = length(nums), ncol = ncol(degs))) 

for(i in 1:length(nums)){
  j <- which(degs[,"Gene.symbol"] == names(nums[i]))
  fac <- degs[j,"logFC"] #(cause we want max LogFC)
  k <- which(max(fac) == fac) #"max" can be changed by the other functions (e.g. mean, min, median and etc.)
  degs.ag[i,] <- degs[j[k],]
}

colnames(degs.ag) <- colnames(degs)
rownames(degs.ag) <- degs.ag[,"Gene.symbol"]

#Find how much upregulate and downregulate genes you had:
degsup <- subset(degs.ag, logFC > 2 & adj.P.Val < 0.05) #OncoGene
degsdown <- subset(degs.ag, logFC < -2 & adj.P.Val < 0.05) #TumorSUpressor

#"degs.ag" is what we need for further analysis:
write.table(degs.ag, file="Degs - aggregated in R.txt", quote=F, sep="\t",
            row.names = F)

#Write down gene symbols for STRING input:
write.table(degs.ag$Gene.symbol, file="Degs GeneNames.txt", quote=F,
            row.names = F, col.names = F)

#Save Data:
setwd("E:/Projects/Mammad Crush CGC/RData")
save(degs,degs.ag,degsdown,degsup,design,tT, file = "DEGs.RData")

#--------------------- Downstream Analysis ---------------------
#Leze bOoOo:


#--------------------- VolcanoPlots for DEGs ---------------------

setwd("E:/Projects/Mammad Crush CGC/plots/")
volc <- subset(tT)

volc$Significant <- "No"
volc$Significant[volc$logFC > 2 & volc$adj.P.Val < 0.05] <- "Up"
volc$Significant[volc$logFC < -2 & volc$adj.P.Val < 0.05] <- "Down"

#tiff("volcano plot.tiff",height = 1600, width = 2000,res=300,units = "px" )
pdf("volcano.pdf")

ggplot(volc, aes(logFC, -log10(adj.P.Val),color =Significant))+ 
  geom_point(size = 1.6,shape = 19) + theme_bw()+ 
  geom_vline(xintercept=c(-2, 2), col="#EF00FF",linetype = "longdash") + 
  geom_hline(yintercept=-log10(0.05), col="#EF00FF",linetype = "longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot")



ggplot(volc, aes(logFC, -log10(adj.P.Val), color = Significant)) + 
  geom_point(size = 1.6, shape = 19) + 
  theme_bw() + 
  geom_vline(xintercept=c(-2, 2), col="#EF00FF",linetype = "longdash") + 
  geom_hline(yintercept=-log10(0.05), col="#EF00FF",linetype = "longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot") +
  geom_text(aes(label = Gene.symbol), size = 1.5, nudge_y = 0.2)


dev.off()


#--------------------- HEATMAP for DEGs ---------------------
#-------------- Heatmapoo is an expression matrix of my DEGs (also for WGCNA)
#--- NEW WAY! ---#
#make a dataframe of expression matrix:
bayan <- data.frame(ex)
# extract expression of DEGs for downstream analysis
heatmapoo <- bayan[rownames(bayan) %in% degs.ag$ID,]

#replace probe IDs in rownames heatmapoo into gene symbol:
# create a named vector of gene symbols from degs.ag
symbol_match <- setNames(degs.ag$Gene.symbol, degs.ag$ID)
# replace rownames in heatmapoo with gene symbols
rownames(heatmapoo) <- symbol_match[rownames(heatmapoo)]

#Save Data:
setwd("E:\\Projects\\CGC2023/CGC Mammad_Esophageal Cancer/plots/")
save(heatmapoo, file = "HEATMAPOO.RData")
#--------------#

#making annotation_col:
coloo <- colnames(heatmapoo) 
annotation1 <- data.frame(gr, coloo)
rownames(annotation1) <- annotation1$coloo
colnames(annotation1) = "Group"
annotation1 <- annotation1[,-2,drop=F]


annotation2 <- data.frame(degs.ag$logFC)
rownames(annotation2) <- degs.ag$Gene.symbol
colnames(annotation2) = "Up/Down regulated"
annotation2$`Up/Down regulated`[degs.ag$logFC > 2 & degs.ag$adj.P.Val < 0.05] <- "Up"
annotation2$`Up/Down regulated`[degs.ag$logFC < -2 & degs.ag$adj.P.Val < 0.05] <- "Down"


ann_colors = list(
  "Group" = c( Paratumor= "#7570B3",  Tumor= "#3CF20F"),
  "Up/Down regulated" = c(Down = "#D7A3E6", Up = "#E7298A"))
#get the color code from rapidtables RGB colors
#https://www.rapidtables.com/web/color/RGB_Color.html

setwd("E:\\Projects\\CGC2023/CGC Mammad_Esophageal Cancer/plots/")

png("heatmap_default_2.png",height = 1600, width = 2000,res=300,units = "px" )

pheatmap(heatmapoo,color= bluered(256), show_rownames=F, show_colnames = F,
         border_color =NA,annotation_colors = ann_colors,
         annotation_col = annotation1, annotation_row = annotation2, 
         scale = "row" , angle_col=45,fontsize_col=4.5, fontsize_row = 4.5,
         annotation_names_row = F,annotation_names_col = F,
         clustering_distance_rows = "correlation", #Pearson correlation
         clustering_distance_cols= "correlation",  #Pearson correlation
         clustering_method = "complete")


dev.off()


#--------------------- CLUSTER PROFILER ---------------------

setwd("E:/Projects/Mammad Crush CGC/ClusterP/")
#remotes::install_github("YuLab-SMU/clusterProfiler") 
#https://www.biostars.org/p/305258/
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DOSE)
library(cowplot)

#Load list of genes as character:
ent_gene <- degs.ag$Gene.ID

########## ClusterProfiler (GO & KEGG) ########## 
## GO enrichment 
ego <- enrichGO(gene = ent_gene,
                OrgDb = org.Hs.eg.db ,
                ont = "ALL" ,
)

#use slice function to select top 10 in each 3 GO terms:(and make a dataframe)

egom <-arrange(ego@result, p.adjust) %>% 
  group_by(ONTOLOGY) %>% 
  slice(1:10)

#cALCULATE RICH FACTOR:
egom <- mutate(egom, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#write enrichGO
ego_table <- mutate(ego@result, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

write.table(ego_table, file="enrichGO.txt", quote=F, sep="\t")

####### dotplot BY Rich factor:

d1 <- ggplot(egom, 
             aes(richFactor,Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Adjusted P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  
  guides(size=guide_legend("Gene number",order = 1),
         "Adjusted P-Value" = guide_legend(order = 2))+
  ggtitle("Statistics of GO Enrichment")+
  facet_grid(rows = vars(ONTOLOGY),scales = "free_y")


## KEGG enrichment :
options(clusterProfiler.download.method = "wininet") # cause KEGG API has moved to HTTPS 
#(https://github.com/YuLab-SMU/clusterProfiler/issues/470)

#gene_vector <- setNames(rep(1, length(ent_gene)), ent_gene)

ekg <- enrichKEGG(ent_gene,
                  organism = "hsa",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)

#deactive pvalue and qvalue effect on the numbers of pathways in results (get all results)

view(ekg)

#view gene symbols:

ekgm<- setReadable(ekg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#cALCULATE RICH FACTOR:
ekgm <- mutate(ekgm, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#make data frame to read by ggplot2:
ekgm <- as.data.frame(ekgm)
#write enrichKEGG
write.table(ekgm, file="enrichKEGG.txt", quote=F, sep="\t")


####### dotplot BY Rich factor:

# select first 20:
ekgm <-arrange(ekgm, p.adjust) %>% slice(1:20)

d2 <- ggplot(ekgm,
             aes(richFactor,Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Adjusted P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  guides(size=guide_legend("Gene number"))+
  ggtitle("Statistics of KEGG Enrichment")


#COMBINE DOTPLOT:
tiff("dotplot_combine.tiff",width = 5500, height = 2300,res = 300,units = "px")

plot_grid(
  d2 ,d1,
  rel_widths = c(1.75, 2),
  nrow = 1,
  labels = c('A', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "bold",
  label_colour = "black",
  label_size = 30
)

dev.off()






