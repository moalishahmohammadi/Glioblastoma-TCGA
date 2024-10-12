setwd("X:\\Projects\\2nd congressTCGA\\Data")

#### Loading required packages
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readxl)
library(DT)
library(edgeR)
library(dplyr)

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
#GBM.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(GBMMatrix)






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

###Enrichment analysis
library(enrichplot)

#Makingg my genelist

)









