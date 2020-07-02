# Size of the figures should be a single (86 mm) or double (178 mm) column width.
library(ggplot2)
library(statsExpressions)

load('FileContent.RData')

# Number of cells
nrow(fContent)
# Number of experiments
length(levels(fContent$SRS))
# Number of projects
length(levels(fContent$SRA))

# Extracting Information
MT <- fContent
table(MT$MTCOUNT != 0)
MT <- MT[MT$MTCOUNT != 0,]
MT <- as.data.frame.array(MT)
MT <- MT[order(MT$TOTALCOUNT),]
MT$log10MTCOUNT <- log10(MT$MTCOUNT)
MT$log10TOTALCOUNT <- log10(MT$TOTALCOUNT)
MT$log10NGENES <- log10(MT$NGENES)

# nGenes regression
expectedGENE <- predict(lm(MT$log10NGENES~poly(MT$log10TOTALCOUNT,2)), newdata = data.frame(TOTALCOUNT = MT$log10TOTALCOUNT), interval = 'prediction', level = 0.95)
expectedGENE <- as.data.frame(expectedGENE)
MT$expGENElwr <- expectedGENE$lwr
MT$expGENE <- expectedGENE$fit
MT$expGENEupr <- expectedGENE$upr

MT$CAT <- paste0(MT$PROTOCOL, ' - ', MT$SP)
CORVAL <- cor(MT$NGENES,MT$TOTALCOUNT, method = 'sp')
cor.test(MT$NGENES,MT$TOTALCOUNT, method = 'sp')

MT <- MT[MT$log10NGENES > MT$expGENElwr & MT$log10NGENES < MT$expGENEupr,]

# Mitochondrial counts regression
expectedMT <- predict(lm(MT$log10MTCOUNT~MT$log10TOTALCOUNT), newdata = data.frame(TOTALCOUNT = MT$log10TOTALCOUNT), interval = 'prediction', level = 0.95)
expectedMT <- as.data.frame(expectedMT)
MT$expMTCOUNTlwr <- expectedMT$lwr
MT$expMTCOUNT <- expectedMT$fit
MT$expMTCOUNTupr <- expectedMT$upr

CORVAL <- cor(MT$log10TOTALCOUNT, MT$log10MTCOUNT)
cor.test(MT$log10TOTALCOUNT, MT$log10MTCOUNT)

table(MT$log10MTCOUNT > MT$expMTCOUNTlwr & MT$log10MTCOUNT < MT$expMTCOUNTupr)
table(MT$log10MTCOUNT > MT$expMTCOUNTlwr)
table(MT$log10MTCOUNT < MT$expMTCOUNTupr)
MT <- MT[MT$log10MTCOUNT > MT$expMTCOUNTlwr & MT$log10MTCOUNT < MT$expMTCOUNTupr,]

# Comparison by specie
p1 <- ggplot(MT, aes(x=SP, y=MTRATIO)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  coord_flip() + geom_jitter(shape=16, position=position_jitter(0.1), cex = 0.05, alpha = 0.01) +
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') +
  #geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Specie') +
  ylab('') + 
  labs(tag = 'A')

# Comparison by technology
p2 <- ggplot(MT, aes(x=CAT, y=MTRATIO)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) +
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') +
  #geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Technology and Specie') +
  ylab('') + 
  labs(tag = 'B')

# Human Tissues
C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Homo sapiens' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 1000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'
C10X$TISSUE[C10X$TISSUE == 'Peripheral blood mononuclear cell'] <- 'PBMC'
C10X$TISSUE[C10X$TISSUE == 'Peripheral blood mononuclear cells'] <- 'PBMC'
tT <- table(C10X$TISSUE)
humanTISSUES <- t(sapply(names(tT), function(X){
  summary(C10X$MTRATIO[C10X$TISSUE == X], )
}))
write.csv(round(humanTISSUES,3), 'mtProportion_HumanTissues.csv')

C10X$TISSUE <- as.factor(C10X$TISSUE)
levels(C10X$TISSUE) <- paste0(names(tT), ' (', tT, ')')

mT <- sapply(unique(C10X$TISSUE), function(X){
  median(C10X$MTRATIO[C10X$TISSUE %in% X])
})
names(mT) <- unique(C10X$TISSUE)
C10X$TISSUE <- factor(C10X$TISSUE, levels = names(sort(mT)))

p3 <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Human Tissue') +
  ylab('Mitochondrial Ratio') + labs(tag = 'C')

# Mouse Tissues
C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Mus musculus' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 1000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'
C10X$TISSUE[C10X$TISSUE == 'Cortex 1'] <- 'Cortex'
C10X$TISSUE[C10X$TISSUE == 'Cortex 2'] <- 'Cortex'
C10X$TISSUE[C10X$TISSUE == 'Cortex 3'] <- 'Cortex'

tT <- table(C10X$TISSUE)
mouseTISSUES <- t(sapply(names(tT), function(X){
  summary(C10X$MTRATIO[C10X$TISSUE == X], )
}))
write.csv(round(mouseTISSUES,3), 'mtProportion_MouseTissues.csv')
C10X$TISSUE <- as.factor(C10X$TISSUE)

levels(C10X$TISSUE) <- paste0(names(tT), ' (', tT, ')')

mT <- sapply(unique(C10X$TISSUE), function(X){
  median(C10X$MTRATIO[C10X$TISSUE %in% X])
})
names(mT) <- unique(C10X$TISSUE)
C10X$TISSUE <- factor(C10X$TISSUE, levels = names(sort(mT)))


p4 <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 10) + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  xlab('Mouse Tissue') +
  ylab('Mitochondrial Ratio') + labs(tag = 'D')

library(patchwork)
png('Figures/M1.png', width = 4000, height = 3500, res = 300)
layout <- 'AABB
CCBB
CCBB
DDBB
DDBB
DDBB
DDBB
DDBB
DDBB
DDBB'
p1 + p4 + p2 + p3 + plot_layout(design = layout, )
dev.off()

