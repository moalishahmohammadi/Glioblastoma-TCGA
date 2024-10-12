setwd("X:\\Projects\\2nd congressTCGA\\Data")

#### Loading required packages
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readxl)
library(DT)
library(edgeR)
library(dplyr)
library(limma)

#### Characterizing and downloading the desired data

query.exp <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts"
                      )

GDCdownload(query = query.exp , method = "api" , files.per.chunk = 5 )


# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
exp <- GDCprepare(query.exp , summarizedExperiment = T, )
GBMMatrix <- assay(exp,"unstranded")
write.csv(GBMMatrix,"GBMMatrix.csv")
#Extracting Clinical info as metadata
clinical<- colData(exp)
#transfering clinical status into column bar

#selecting specific columns
clinical<- as.data.frame(clinical[ ,c(1, 2, 8, 10, 12, 13, 17, 22, 29, 36, 41, 43)])
write.table(clinical, "clinical.txt")
###TCGAanalyze_DEA & TCGAanalyze_LevelTab: Differential expression analysis (DEA)###
edgeR::DGEList 
edgeR::estimateCommonDisp
edgeR::exactTest 
edgeR::topTags
#This function receives as arguments:
#mat2 The matrix of the second group (in the example, group 2 is tumor samples)
#Cond1type Label for group 1
#Cond1type Label for group 2

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = GBMMatrix, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 2,
  method = "glmLRT"
)

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)

write.csv(dataDEGsFiltLevel,"dataDEGsFiltLevel.csv")

# Load the biomaRt package
library(biomaRt)

# Define the Ensembl dataset you want to use (e.g., human genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get the mapping for gene IDs to gene names
geneIDToGeneNameMap <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)

# Load your DEGs dataframe (dataDEGsFiltLevel)
# Replace "your_dataDEGsFiltLevel.csv" with your actual data file path
dataDEGsFiltLevel <- read.csv("dataDEGsFiltLevel.csv")

# Remove the 'X' column
dataDEGsFiltLevel <- dataDEGsFiltLevel[, !names(dataDEGsFiltLevel) %in% "X"]

# Merge the mapping with your DEGs dataframe using the 'mRNA' column
dataDEGsFiltLevel <- merge(dataDEGsFiltLevel, geneIDToGeneNameMap, by.x = "mRNA", by.y = "ensembl_gene_id", all.x = TRUE)

# Rename the 'external_gene_name' column to 'X' (optional)
colnames(dataDEGsFiltLevel)[colnames(dataDEGsFiltLevel) == "external_gene_name"] <- "X"

# Remove the original 'mRNA' column
dataDEGsFiltLevel <- dataDEGsFiltLevel[, !names(dataDEGsFiltLevel) %in% "mRNA"]

# Save the updated dataframe to a CSV file
write.csv(dataDEGsFiltLevel, file = "updated_dataDEGsFiltLevel.csv", row.names = FALSE)

# Rename the 'X' column to 'name'
colnames(dataDEGsFiltLevel)[colnames(dataDEGsFiltLevel) == "X"] <- "Gene"

# Move the 'name' column to the first position
dataDEGsFiltLevel <- dataDEGsFiltLevel[, c("Gene", setdiff(colnames(dataDEGsFiltLevel), "Gene"))]



######################################################
# Assuming "GeneName" is the column with gene names in your geneIDToGeneNameMap dataframe
duplicate_genes <- geneIDToGeneNameMap$GeneName[duplicated(geneIDToGeneNameMap$GeneName)]

# Check if there are any duplicates
if (length(duplicate_genes) > 0) {
  cat("Duplicate gene names found:\n")
  print(duplicate_genes)
} else {
  cat("No duplicate gene names found.\n")
}


######################################################

setwd("X:\\Projects\\2nd congressTCGA\\Plots")
#Let's Visualize
#Volcano Plot
TCGAVisualize_volcano(x = dataDEGsFiltLevel$logFC,
                      y = dataDEGsFiltLevel$FDR,
                      x.cut = 2,
                      y.cut = 0.05,
                      width = 15,
                      height = 10,
                      legend = "State",
                      color = c("black","blue","pink"),
                      xlab = "Gene expression fold change (Log2)",
                      title = "Volcano plot (Primary solid Tumor vs Solid Tissue Normal)",
                      filename = "Volcano plot of STAD-DEGs.pdf",
                      show.names = F) 
dev.off()

library(TCGAbiolinks)

Genelist <- listam

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)

setwd("X:\\Projects\\2nd congressTCGA\\Plots")


TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP), 
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = Genelist, 
  nBar = 10
)










library(TCGAbiolinks)
library(ggplot2)

# Define your list of genes (replace listam with your actual list)
Genelist <- listam  # Replace listam with your actual list of genes

# Perform Enrichment Analysis
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)

# Extract the results for visualization
enrichedGO <- ansEA$ResBP  # Replace with the appropriate results table
enrichedPathways <- ansEA$ResPat  # Replace with the appropriate results table

# Visualize the enriched GO terms using ggplot2
ggplot(enrichedGO, aes(x = Description, y = GeneRatio, fill = pvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Term", y = "Gene Ratio", fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Visualize the enriched pathways using ggplot2
ggplot(enrichedPathways, aes(x = Pathway, y = GeneRatio, fill = pvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Pathway", y = "Gene Ratio", fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))













###Bar Plot

library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

library(enrichplot)
barplot(edo, showCategory=20) 

dev.off()


###Dot Plot

png("dotplot.png", width = 800, height = 900, res = 100)

edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

dev.off()

###gene concept network

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 



pdf("cowplot.pdf")
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
dev.off()


p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 

pdf("cowplot2.pdf")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()



###tree plot

pdf("treeplot.pdf")
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')
dev.off()


