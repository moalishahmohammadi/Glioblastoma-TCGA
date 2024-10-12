setwd("X:\\Projects\\2nd congressTCGA\\Data")

### removal of control samples

RNAseq = GBMMatrix[apply(GBMMatrix,1,function(x) sum(x==0))<ncol(GBMMatrix)*0.8,]


library(limma)
RNAseq_voom = voom(RNAseq)$E

#excluding junkgenes or genes with low counts
#transpose matrix to correlate genes in the following
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])


#finiding out nodes which have lots of connection
#similarity measure between gene profiles: biweight midcorrelation
library(WGCNA)
s = abs(bicor(WGCNA_matrix))


# picking the right power to generate network
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')


#calculation of adjacency matrix
beta = 3
a = s^beta

#dissimilarity measure
w = 1-a


#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)

#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')







library(ape)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))











#load clinical metadata. Make sure that patient barcodes are in the same format 
#create second expression matrix for which the detailed clinical data is available 
WGCNA_matrix2 = WGCNA_matrix[match(clinical$Name, rownames(WGCNA_matrix)),]

#CAVE: 1 sample of detailed clinical metadata is not in downloaded data (TCGA-GN-A269-01')
not.available = which(is.na(rownames(WGCNA_matrix2))==TRUE)
WGCNA_matrix2 = WGCNA_matrix2[-not.available,]
str(WGCNA_matrix2)

#hence it needs to be removed from clinical table for further analysis
clinical = clinical[-not.available,]









#grouping in high and low lymphocyte score (lscore)
lscore = as.numeric(clinical$LYMPHOCYTE.SCORE)
lscore[lscore<3] = 0
lscore[lscore>0] = 1

#calculate gene significance measure for lymphocyte score (lscore) - Welch's t-Test
GS_lscore = t(sapply(1:ncol(WGCNA_matrix2),function(x)c(t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$p.value,
                                                        t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$estimate[1],
                                                        t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$estimate[2])))
GS_lscore = cbind(GS.lscore, abs(GS_lscore[,2] - GS_lscore[,3]))
colnames(GS_lscore) = c('p_value','mean_high_lscore','mean_low_lscore',
                        'effect_size(high-low score)'); rownames(GS_lscore) = colnames(WGCNA_matrix2)