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

sID <- 130
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

downloadedCells <- subset(downloadedCells, idents = 'Alveolar macrophages')
Idents(downloadedCells) <- downloadedCells$panglaoCluster
B <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Alveolar Macrophages Clusters', subtitle = sampleList$SRS[sID], tag = 'B')

mtCounts <- downloadedCells@assays$RNA@counts[grepl('MT-',rownames(downloadedCells@assays$RNA@counts), ignore.case = TRUE),]
mtCounts <- colSums(mtCounts)
mtProportion <- mtCounts/colSums(downloadedCells@assays$RNA@counts)
downloadedCells$mtProportion <- mtProportion

dF <- data.frame(C = downloadedCells$panglaoCluster, MT = downloadedCells$mtProportion)

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
  labs(title = 'Mitochondrial Proportion', subtitle = paste0('Alveolar macrophages ',sampleList$SRS[sID]), tag = 'C')


Idents(downloadedCells) <- downloadedCells$panglaoCluster
D <- FeaturePlot(downloadedCells, 'mtProportion', reduction = 'tsne', order = TRUE, label = TRUE, repel = TRUE) +
  theme_bw() +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Alveolar Macrophages Mitochondrial Proportion', subtitle = sampleList$SRS[sID], tag = 'D') +
  theme(plot.title = element_text(face = 2))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '1', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
E <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 1 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'E')

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '0', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
F <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 0 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '3', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
G <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 3 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '2', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
H <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 2 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))


png('Figures/alveolarMacrophages.png', width = 4000, height = 2800, res = 300)
(A + B + C) / (D | (E + F +G + H))
dev.off()
