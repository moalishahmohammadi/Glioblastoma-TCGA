setwd("G:\\OOF\\2nd Congress\\R codes")

#### Loading required packages
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readxl)
library(DT)

#### Characterizing and downloading the desired data

query.exp <- GDCquery(project = "TCGA-CESC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",
                      legacy = FALSE)

query.met <- GDCquery(project = "TCGA-CESC",
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.format = "TXT",
                      legacy = FALSE)

GDCdownload(query = query.exp , method = "api" , files.per.chunk = 5 )
GDCdownload(query = query.met , method = "api" , files.per.chunk = 5 )

exp <- GDCprepare(query.exp , summarizedExperiment = T, )

met <- GDCprepare(query.met , summarizedExperiment = T, )


#################################################
####                                         ####
####     Gene Expression Data Analysis       ####
####                                         ####
#################################################

#### Determining tumoral and non-tumoral samples

exp.N <- colData(exp)$definition %in% c("Solid Tissue Normal")
exp.T <- colData(exp)$definition %in% c("Primary solid Tumor")


#### Pre-processing and normalizing expression data

exp.prep <- TCGAanalyze_Preprocessing(object = exp, cor.cut = 0.6)

exp.norm <- TCGAanalyze_Normalization(tabDF = exp.prep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

exp.filt <- TCGAanalyze_Filtering(tabDF = exp.norm,
                                  method = "quantile",
                                  qnt.cut = 0.25)

#### Extraction of differentially-expressed genes

exp.DEGs <- TCGAanalyze_DEA(mat1 = exp.filt[,exp.T],
                            mat2 = exp.filt[,exp.N],
                            Cond1type = "Primary solid Tumor",
                            Cond2type = "Solid Tissue Normal",
                            pipeline = "edgeR",
                            method = "exactTest",
                            fdr.cut = 1,
                            logFC.cut = 0)

write.table(exp.DEGs, "Differentially-expressed Genes.txt", quote = F)

#### Drawing the volcano plot

TCGAVisualize_volcano(x = exp.DEGs$logFC,
                      y = exp.DEGs$FDR,
                      x.cut = 1,
                      y.cut = 0.05,
                      width = 15,
                      height = 10,
                      legend = "State",
                      color = c("black","red","green"),
                      xlab = "Gene expression fold change (Log2)",
                      title = "Volcano plot (Primary solid Tumor vs Solid Tissue Normal)",
                      filename = "Volcano plot of STAD-DEGs.pdf",
                      show.names = F) 



#################################################
####                                         ####
####        Methylation data analysis        ####
####                                         ####
#################################################

#### Pre-processing 
#### Removing NAs, sex chromosomes and SNPs

met <- subset(met, subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met, subset = !as.character(seqnames(met))
              %in% c("*", "chrX", "chrY"))

#### Extraction of differentially-methylated CpGs.

met.state <- TCGAanalyze_DMC(met,
                             groupCol = "definition",
                             group1 = "Primary solid Tumor",
                             group2 = "Solid Tissue Normal",
                             p.cut = 0.05,
                             diffmean.cut = 0.2,
                             legend = "State",
                             plot.filename = "Volcano plot of STAD-DMCs.pdf")

met.DMCs <- subset(met.state,
                   subset = met.state$status != "Not Significant")

met.DMC.hypo <- subset(met.state,
                       met.state$status == "Hypomethylated in Primary solid Tumor")

met.DMC.hyper <- subset(met.state,
                        met.state$status == "Hypermethylated in Primary solid Tumor")

#################################################
####                                         ####
####      Annotation of expression data      ####        
####                                         ####
#################################################

genes<- read.delim("gene-report.txt")

genes$<-substr(genes$gene_id, 1 , 15 )

genes=aggregate(x =genes,by=list(genes$x), max )

row.names(genes)<- genes$x


exp.DEGs$gene_name<- genes[as.character(rownames(exp.DEGs)),"gene_name"]

exp.DEGs$gene_type<- genes[as.character(rownames(exp.DEGs)),"gene_type"]

write.table(exp.DEGs,"diff with genes name.txt",
            quote = F,sep = "\t")




#################################################
####                                         ####
####      Annotation of methylation data     ####      
####                                         ####
#################################################

annot.met <- as.data.frame(read.delim("annotation.txt"))

annot.met<-annot.met[,c("Composite.Element.REF","Gene_Symbol")]

dmc<- as.data.frame(met.state[,c("mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal",
                                 "p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal")])

colnames(dmc)<-c("delta beta","adj pvalue")



data_common <- merge(dmc,annot.met,by.x="row.names",by.y="Composite.Element.REF")


data_common <- tidyr::separate_rows(data_common, Row.names, Gene_Symbol,sep = ";") 

data_common$comb <- paste(data_common$Row.names,data_common$Gene_Symbol)

data_common <- data_common[!duplicated(data_common$comb), ]

data_common <- data_common[, -5]

write.table(data_common,"data_common.txt",
            quote = F,sep = "\t")
#################################################
####                                         ####
####      integration and visualization      ####           
####                                         ####
#################################################
#adj.p,symbol,logFc,deltaBeta,symbol

deg<- read.delim("F:/TCGA_yasaman/diff with genes name.txt")


dmg<- read.delim("F:/TCGA_yasaman/data_common.txt")


deg<-deg[,c(2,5,6)]

dmg<-dmg[,c(2,3,4)]

dmg1=aggregate(x =dmg,by=list(dmg$Gene_Symbol), max )
deg1=aggregate(x =deg,by=list(deg$Gene_Symbol), max )


# names in dmg are different . they have more space so delete the space by below cod

dmg$Gene_Symbol <- gsub(' ','',dmg$Gene_Symbol)

data<-merge(deg,dmg,by.x= "Gene_Symbol" , by.y= "Gene_Symbol")


data$Group <- ifelse(data$logFC <(-1) & data$`delta.beta` <(-0.2) & 
                       data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hypo-Down" ,
                     ifelse(data$logFC >1 & data$`delta.beta` <(-0.2) &
                              data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hypo-Up",
                            ifelse(data$logFC <(-1) & data$`delta.beta` >0.2 &
                                     data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hyper-Down",
                                   ifelse(data$logFC >1 & data$`delta.beta` >0.2 &
                                            data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hyper-Up",
                                          "Not-sig"))))


#write.csv(data,"final.csv")
#x<-data[data$Group !="Not-sig",]

write.table(data,"data.txt",
            quote = F,sep = "\t")
# plotting

library(ggplot2)

cols <- c("Hypo-Down" = "#B8860B", "Hypo-Up" = "blue", "Not-sig" = "grey", "Hyper-Down" = "red", "Hyper-Up" = "springgreen4")

ggplot(data, aes(x = `delta.beta`, y =logFC , color = Group)) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = cols) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_hline(yintercept = -1, linetype="dashed") + 
  geom_vline(xintercept = 0.2, linetype="dashed") + 
  geom_vline(xintercept = -0.2,  linetype="dashed") +
  xlab("Mean methylation diffrences") + 
  ylab("Log2 expression change") 


#################################################
####                                         ####
####      cluster profiler                   ####           
####                                         ####
#################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DOSE)

## GO enrichment 
### dot GO:

setwd("F:/TCGA_yasaman/")

y <- read.delim("F:/TCGA_yasaman/for cluster.txt", sep="\t",header = T)

y2 <- y$gene_name
ent_uni <- bitr(y2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ent_uni <- ent_uni$ENTREZID
ego_BP <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "BP" ,
                   #universe = ent_uni 
)

ego_MF <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "MF" ,
                   #universe = ent_uni 
)

ego_CC <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "CC" ,
                   #universe = ent_uni 
)



#N <- as.numeric(sub("\\d+/", "", ego[1, "BgRatio"]))


p1 <- dotplot(ego_BP, showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Biological Process") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))



p2 <- dotplot(ego_MF, showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Molecular Function") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))


p3 <- dotplot(ego_CC,  showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Cellular Component") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))


png("dotplot_yasaman.png",width = 2500, height = 4000,res = 300,units = "px")


cowplot::plot_grid(p1, p2, p3, nrow =3,
                   labels=c("I","II","III"))

dev.off()


### cnet GO:

geneList <- y$logFC
names(geneList) <- as.character(y$symbol)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)


edox1 <- setReadable(ego_BP, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)

edox2 <- setReadable(ego_MF, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)

edox3 <- setReadable(ego_CC, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)


k1 <- cnetplot(edox1,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Biological Process") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


k2 <- cnetplot(edox2,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Molecular Function") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


k3 <- cnetplot(edox3,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Cellular Component") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


png("cnetplot_GO.png",width = 10000, height = 3000,res = 300,units = "px")

cowplot::plot_grid(k1, k2, k3, ncol =3,
                   labels=c("I","II","III"))

dev.off()

####################################################



## KEGG enrichment :

ekg <- enrichKEGG(ent_uni,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

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

#tiff("dotplot_KEGG.tiff",width = 3000, height = 2000,res = 300,units = "px")

d2 <- ggplot(ekgm,
             aes(richFactor,Description)) + 
  #geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE),name="Corrected P-Value") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  #theme_minimal() + 
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  guides(size=guide_legend("Gene number"))+
  ggtitle("Statistics of KEGG Enrichment")

#dev.off()

####################################################

### cnet KEGG:
geneList <- x$logFC
names(geneList) <- as.character(x$geneSymbol)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)


jj <-arrange(ekg@result, p.adjust)
ekg@result <- jj

ekgx <- setReadable(ekg, 'org.Hs.eg.db', 'ENTREZID')
view(ekgx)


#tiff("cnetplot_KEGG.tiff",width = 4300, height = 3500,res = 300,units = "px")

p2 <- cnetplot(ekgx,colorEdge = TRUE,categorySize=2,
               foldChange = geneList,circular = TRUE) + 
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold"))


#dev.off()

#COMBINE CNET:
tiff("cnetplot_combine.tiff",width = 8500, height = 3500,res = 300,units = "px")

plot_grid(
  p2 ,NULL,p1,
  rel_widths = c(1, 0.05, 1),
  nrow = 1,
  labels = c('A','', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_colour = "black",
  label_size = 30
)

dev.off()

#COMBINE DOTPLOT:
tiff("dotplot_combine.tiff",width = 5500, height = 2300,res = 300,units = "px")

plot_grid(
  d2 ,d1,
  rel_widths = c(1.75, 2),
  nrow = 1,
  labels = c('A', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_colour = "black",
  label_size = 30
)

dev.off()
##############################

### Enrichment Map using KEGG :
ekg <- enrichKEGG(ent_gene,
                  pvalueCutoff = 1)



tiff("EnrichmentMap_KEGG.tiff",width = 4000, height = 4000,res = 300,units = "px")

emapplot(ekg,showCategory = 30,pie="count", pie_scale=1.5, layout="kk")

dev.off()



