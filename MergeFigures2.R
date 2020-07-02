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
  labs(title = 'Cardiomyocytes Mitochondrial Proportion', subtitle = sampleList$SRS[sID], tag = 'B') +
  theme(plot.title = element_text(face = 2)) + 
  xlim(c(-45,15)) +
  ylim(c(-45,15))

DE <- FindMarkers(downloadedCells, ident.1 = '19', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
E <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('19 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'D')

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


png('Figures/Cardiomyocytes.png', width = 6000, height = 1500, res = 300)
layout <- '
AAABBBCCCDDEE
AAABBBCCCFFGG
'
A + D + C + E + F + G + H + plot_layout(design = layout)
dev.off()