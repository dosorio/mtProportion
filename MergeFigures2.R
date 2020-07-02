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
A1 <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cell types', subtitle = sampleList$SRS[sID], tag = 'A1')

downloadedCells <- subset(downloadedCells, idents = 'Cardiomyocytes')
Idents(downloadedCells) <- downloadedCells$panglaoCluster

mtCounts <- downloadedCells@assays$RNA@counts[grepl('MT-',rownames(downloadedCells@assays$RNA@counts), ignore.case = TRUE),]
mtCounts <- colSums(mtCounts)
mtProportion <- mtCounts/colSums(downloadedCells@assays$RNA@counts)
downloadedCells$mtProportion <- mtProportion

dF <- data.frame(C = downloadedCells$panglaoCluster, MT = downloadedCells$mtProportion)
# dF <- dF[downloadedCells$CellTypes %in% 'Cardiomyocytes',]

cID <- sort(unique(dF$C))
mtMedian <- sapply(cID, function(X){median(dF$MT[dF$C %in% X])})

dF$C <- factor(dF$C, levels= cID[order(mtMedian, decreasing = TRUE)])
C1 <- ggplot(dF, aes(MT, C)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(cex = 0.1, alpha = 0.25) +
  theme_bw() +
  geom_vline(xintercept = 0.05, col='red', lty = 2) +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  ylab('Cluster') +
  xlab('Mitochondrial Proportion') +
  labs(title = sampleList$SRS[sID]) +
  labs(title = 'Cardiomyocytes', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'C1')


Idents(downloadedCells) <- downloadedCells$panglaoCluster
D1 <- FeaturePlot(downloadedCells, 'mtProportion', reduction = 'tsne', order = TRUE, label = TRUE, repel = TRUE) +
  theme_bw() +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cardiomyocytes', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'B1') +
  theme(plot.title = element_text(face = 2)) + 
  xlim(c(-45,15)) +
  ylim(c(-45,15))

DE <- FindMarkers(downloadedCells, ident.1 = '19', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
E1 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('19 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'D1')

DE <- FindMarkers(downloadedCells, ident.1 = '0', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
F1 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('0 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '9', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
G1 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('9 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '7', ident.2 = '4', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
H1 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('7 vs 4 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))


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
A2 <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cell types', subtitle = sampleList$SRS[sID], tag = 'A2')

downloadedCells <- subset(downloadedCells, idents = 'Alveolar macrophages')
Idents(downloadedCells) <- downloadedCells$panglaoCluster

mtCounts <- downloadedCells@assays$RNA@counts[grepl('MT-',rownames(downloadedCells@assays$RNA@counts), ignore.case = TRUE),]
mtCounts <- colSums(mtCounts)
mtProportion <- mtCounts/colSums(downloadedCells@assays$RNA@counts)
downloadedCells$mtProportion <- mtProportion

dF <- data.frame(C = downloadedCells$panglaoCluster, MT = downloadedCells$mtProportion)

cID <- sort(unique(dF$C))
mtMedian <- sapply(cID, function(X){median(dF$MT[dF$C %in% X])})

dF$C <- factor(dF$C, levels= cID[order(mtMedian, decreasing = TRUE)])
C2 <- ggplot(dF, aes(MT, C)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(cex = 0.1, alpha = 0.25) +
  theme_bw() +
  geom_vline(xintercept = 0.05, col='red', lty = 2) +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  ylab('Cluster') +
  xlab('Mitochondrial Proportion') +
  labs(title = sampleList$SRS[sID]) +
  labs(title = 'Alveolar macrophages', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'C2')


Idents(downloadedCells) <- downloadedCells$panglaoCluster
D2 <- FeaturePlot(downloadedCells, 'mtProportion', reduction = 'tsne', order = TRUE, label = TRUE, repel = TRUE) +
  theme_bw() +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Alveolar Macrophages', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'B2') +
  theme(plot.title = element_text(face = 2))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '1', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
E2 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 1 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'D2')

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '0', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
F2 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 0 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '3', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
G2 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 3 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '5', ident.2 = '2', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
H2 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('5 vs 2 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

sID <- 56
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
A3 <- TSNEPlot(downloadedCells, label = TRUE, repel = TRUE) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Cell types', subtitle = sampleList$SRS[sID], tag = 'A3')

downloadedCells <- subset(downloadedCells, idents = 'Endothelial cells')
Idents(downloadedCells) <- downloadedCells$panglaoCluster

mtCounts <- downloadedCells@assays$RNA@counts[grepl('MT-',rownames(downloadedCells@assays$RNA@counts), ignore.case = TRUE),]
mtCounts <- colSums(mtCounts)
mtProportion <- mtCounts/colSums(downloadedCells@assays$RNA@counts)
downloadedCells$mtProportion <- mtProportion

dF <- data.frame(C = downloadedCells$panglaoCluster, MT = downloadedCells$mtProportion)

cID <- sort(unique(dF$C))
mtMedian <- sapply(cID, function(X){median(dF$MT[dF$C %in% X])})

dF$C <- factor(dF$C, levels= cID[order(mtMedian, decreasing = TRUE)])
C3 <- ggplot(dF, aes(MT, C)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(cex = 0.1, alpha = 0.25) +
  theme_bw() +
  geom_vline(xintercept = 0.05, col='red', lty = 2) +
  theme(legend.position="none", plot.title = element_text(face = 2)) +
  ylab('Cluster') +
  xlab('Mitochondrial Proportion') +
  labs(title = sampleList$SRS[sID]) +
  labs(title = 'Endothelial Cells', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'C3')


Idents(downloadedCells) <- downloadedCells$panglaoCluster
D3 <- FeaturePlot(downloadedCells, 'mtProportion', reduction = 'tsne', order = TRUE, label = TRUE, repel = TRUE) +
  theme_bw() +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  labs(title = 'Endothelial Cells', subtitle = paste0('Mitochondrial Proportion ', sampleList$SRS[sID]), tag = 'B3') +
  theme(plot.title = element_text(face = 2))

DE <- FindMarkers(downloadedCells, ident.1 = '11', ident.2 = '2', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
E3 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('11 vs 2 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)), tag = 'D3')

DE <- FindMarkers(downloadedCells, ident.1 = '11', ident.2 = '0', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
F3 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('11 vs 0 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '11', ident.2 = '3', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
G3 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('11 vs 3 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))

DE <- FindMarkers(downloadedCells, ident.1 = '11', ident.2 = '7', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
PValue <- fgseaMultilevel(KEGG['Apoptosis'], FC)
H3 <- plotEnrichment(KEGG$Apoptosis, FC) +
  theme_bw() +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'Apoptosis', subtitle = paste0('11 vs 7 | NES = ',round(PValue$NES,2), ' | P = ', formatC(PValue$padj, format = 'e', digits = 2)))



png('Figures/M2.png', width = 6500, height = 1300*3, res = 300)
layout <- '
AAABBBCCCDDEE
AAABBBCCCFFGG
HHHIIIJJJKKLL
HHHIIIJJJMMNN
OOOPPPQQQRRSS
OOOPPPQQQTTUU
'
A1 + D1 + C1 + E1 + F1 + G1 + H1 + A2 + D2 + C2 + E2 + F2 + G2 + H2 + A3 + D3 + C3 + E3 + F3 + G3 + H3 + plot_layout(design = layout)
dev.off()

