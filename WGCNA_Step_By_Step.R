#Co-Expression Analysis using WGCNA

#Packages installation
#BiocManager::install("WGCNA")

#Set Working directory
setwd("X:\\Projects\\2nd congressTCGA\\Data"); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#it can be RNA seq data or micro array
ExpressionData = read.csv("GBMMatrix.csv", row.names = 1)
transposed_ExpressionData = as.data.frame(t(ExpressionData));

#filter high NA values and zero variance genes and samples
gsg = goodSamplesGenes(transposed_ExpressionData, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(transposed_ExpressionData)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(transposed_ExpressionData)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  transposed_ExpressionData = transposed_ExpressionData[gsg$goodSamples, gsg$goodGenes]
}


# Outlier detection using hierarchical clustering
sampleTree = hclust(dist(transposed_ExpressionData), method = "average");
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

{plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# adjust abline to make a better cut
  abline(h = 35, col = "red") }

# Remove outliers here
# Plot a line to show the cut
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 35, minSize = 10)
table(clust)
# select which cluster samples you want to keep.
keepSamples = (clust == 0)
transposed_ExpressionData = transposed_ExpressionData[keepSamples, ]
nGenes = ncol(transposed_ExpressionData)
nSamples = nrow(transposed_ExpressionData)

# Load pheno data
traitData = read.csv("Pheno.csv", row.names = 1);

# to be in same order like transposed expression data
traitData = traitData[rownames(transposed_ExpressionData),]

collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(transposed_ExpressionData), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")

save(transposed_ExpressionData, traitData, file = "dataInput.RData")

#=====================================================================================
#
#  Let's make network
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(from = 30, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(transposed_ExpressionData, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
{plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")}
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#put first value which cut the 0.9  as soft power
softPower = 48;
adjacency = adjacency(transposed_ExpressionData, power = softPower);


#=====================================================================================
#
#  Similarity hierarchical clustering of genes
#
#=====================================================================================
# To minimize effects of noise and spurious associations, 
# we transform the adjacency into Topological Overlap Matrix (Similarity)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=====================================================================================
#
#  Let's get modules
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
#minimum genes in a module
minModuleSize = 15;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Let's merge modules using eigengenes
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(transposed_ExpressionData, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")



# Altitude to merge
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(transposed_ExpressionData, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#Plot merged modules
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
table(moduleLabels)
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep.RData")