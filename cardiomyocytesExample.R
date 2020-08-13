library(Seurat)
library(rPanglaoDB)
library(Matrix)
library(ggplot2)
library(patchwork)
library(fgsea)
library(patchwork)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')


sampleList <- getSampleList()
sampleList$Cells <- as.numeric(sampleList$Cells)
sampleList <- sampleList[order(sampleList$Cells, decreasing = TRUE),]
rownames(sampleList) <- NULL

sID <- 47
downloadedCells <- getSamples(srs = sampleList$SRS[sID])
downloadedCells <- mergeExperiments(downloadedCells)
downloadedCells <- NormalizeData(downloadedCells)
downloadedCells <- ScaleData(downloadedCells)
downloadedCells <- FindVariableFeatures(downloadedCells)
vFeaturesNoMito <- VariableFeatures(downloadedCells)
vFeaturesNoMito <- vFeaturesNoMito[!grepl('Mt-', vFeaturesNoMito, ignore.case = TRUE)]
downloadedCells <- RunPCA(downloadedCells, npcs = 20, features = vFeaturesNoMito, verbose = FALSE)
downloadedCells <- RunTSNE(downloadedCells, dims = 1:20)
downloadedCells$CellTypes[!downloadedCells$CellTypes %in% names(table(downloadedCells$CellTypes)[table(downloadedCells$CellTypes) > 100])] <- NA

Idents(downloadedCells) <- downloadedCells$CellTypes
A <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cell types', subtitle = sampleList$SRS[sID], tag = 'A')

downloadedCells <- subset(downloadedCells, idents = 'Cardiomyocytes')
Idents(downloadedCells) <- downloadedCells$panglaoCluster
B <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cardiomyocytes Clusters', subtitle = sampleList$SRS[sID], tag = 'B') + 
  xlim(c(-45,15)) +
  ylim(c(-45,15))

mtCounts <- downloadedCells@assays$RNA@counts[grepl('MT-',rownames(downloadedCells@assays$RNA@counts), ignore.case = TRUE),]
mtCounts <- colSums(mtCounts)
mtProportion <- mtCounts/colSums(downloadedCells@assays$RNA@counts)
downloadedCells$mtProportion <- mtProportion

dF <- data.frame(C = downloadedCells$panglaoCluster, MT = downloadedCells$mtProportion)
# dF <- dF[downloadedCells$CellTypes %in% 'Cardiomyocytes',]

cID <- sort(unique(dF$C))
mtMedian <- sapply(cID, function(X){median(dF$MT[dF$C %in% X])})

dF$C <- factor(dF$C, levels= cID[order(mtMedian, decreasing = TRUE)])
C <- ggplot(dF, aes(MT, C)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(cex = 0.1, alpha = 0.25) +
  theme_bw() +
  geom_vline(xintercept = 0.05, col='red', lty = 2) +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  ylab('Cluster') +
  xlab('Mitochondrial Proportion') +
  labs(title = sampleList$SRS[sID]) +
  labs(title = 'Cardiomyocytes Mitochondrial Proportion', subtitle = paste0(sampleList$SRS[sID]), tag = 'C')


Idents(downloadedCells) <- downloadedCells$panglaoCluster
D <- FeaturePlot(downloadedCells, 'mtProportion', reduction = 'tsne', order = TRUE, label = TRUE, repel = TRUE) +
  theme_bw() +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cardiomyocytes Mitochondrial Proportion', subtitle = sampleList$SRS[sID], tag = 'D') +
  theme(plot.title = element_text(face = 2)) + 
  xlim(c(-45,15)) +
  ylim(c(-45,15))

DE <- FindMarkers(downloadedCells, ident.1 = '19', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))


#### Power analysis #####
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
set.seed(1)
powerAnalysis <- t(pbapply::pbsapply(seq_len(1000), function(X){
  sapply(seq(0,6057, 500)[-1], function(B){
    as.numeric(fgseaMultilevel(KEGG['Apoptosis'], sample(FC, B))$pval)
  })
}))
powerAnalysis <- apply(powerAnalysis,2,as.numeric)
nGenes <- seq(0,6057, 500)[-1]
meanP <- -log10(colMeans(powerAnalysis, na.rm = TRUE))
sdP <- -log10(apply(powerAnalysis,2,function(X){sd(X, na.rm = TRUE)}))
lbP <- meanP - sdP
lbP[lbP < 0] <- 0
ubP <- meanP + sdP
DF <- data.frame(genes=nGenes, mean = meanP, ub = ubP, lb = lbP)
colnames(powerAnalysis) <- nGenes
paDF <- reshape2::melt(powerAnalysis)
paDF <- paDF[,2:3]
paDF$value <- -log10(paDF$value)
colnames(paDF) <- c('nGenes', 'P')
paDF$nGenes <- as.factor(paDF$nGenes)
png('Figures/Power.png', width = 1000, height = 750, res = 300)
# ggplot(DF, aes(genes,mean)) + ylab(expression(-log[10]~P-value)) + xlab('Number of Tested Genes') +
#   ylim(c(0,5)) + 
#   geom_point() +
#   geom_errorbar(mapping = aes(ymin=lb, ymax=ub)) + 
#   geom_hline(yintercept = -log10(0.05), lty = 2, col = 'red') +
#   labs(title = 'Apoptosis') +
#   theme(plot.title = element_text(face = 2)) +
#   theme_bw()
ggplot(paDF, mapping = aes(x = nGenes, y = P)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.25, pch = 16, cex = 0.2) + theme_bw() +
  labs(title = 'Apoptosis') +
  theme(plot.title = element_text(face = 2)) + 
  geom_hline(yintercept = -log10(0.05), lty = 2, col = 'red') + ylab(expression(-log[10]~Enrichment~P-value)) + xlab('Number of Tested Genes') + coord_flip()
dev.off()



E <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('19 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'E')

DE <- FindMarkers(downloadedCells, ident.1 = '0', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
F <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('0 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '9', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
G <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('9 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '7', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
H <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('7 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))


png('Figures/Cardiomyocytes.png', width = 4000, height = 2800, res = 300)
(A | B | C)/(D | (E + F + G + H))
dev.off()
