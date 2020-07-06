# Size of the figures should be a single (86 mm) or double (178 mm) column width.
library(ggplot2)
library(statsExpressions)

load('FileContent.RData')

# Filtering 10X
fContent <- fContent[fContent$PROTOCOL %in% '10x chromium',]

# Shared Tissues
hsaT <- unique(fContent$TISSUE[fContent$SP %in% 'Homo sapiens'])
mmuT <- unique(fContent$TISSUE[fContent$SP %in% 'Mus musculus'])
sharedT <- intersect(hsaT,mmuT)
sharedT <- sharedT[!sharedT %in% 'Unknown']

# Testing
outT <- sapply(seq_along(sharedT), function(X){
  hsaData <- fContent$MTRATIO[fContent$TISSUE %in% sharedT[X] & fContent$SP %in% 'Homo sapiens']
  mmuData <- fContent$MTRATIO[fContent$TISSUE %in% sharedT[X] & fContent$SP %in% 'Mus musculus']
  if(length(hsaData) >= 1000 & length(mmuData) >= 1000){
    testOut <- t.test(hsaData,mmuData, conf.level = 0.95)
    c(testOut$estimate, testOut$conf.int, testOut$p.value)  
  }
})
outT <- t(outT)
colnames(outT) <- c('mean Hsa', 'mean Mmu', 'lower CI', 'upper CI', 'p-value')           
outT <- as.data.frame(outT)
outT$`p-adj` <- p.adjust(outT$`p-value`, method = 'fdr')
outT$`Hsa > Mmu` <- outT$`mean Hsa` > outT$`mean Mmu`
rownames(outT) <- sharedT

write.csv(outT, 'spCT.csv')


# Shared Cell Types
hsaT <- unique(fContent$CT[fContent$SP %in% 'Homo sapiens'])
mmuT <- unique(fContent$CT[fContent$SP %in% 'Mus musculus'])
sharedCT <- intersect(hsaT,mmuT)
sharedCT <- sharedCT[!sharedCT %in% 'Unknown']


outCT <- sapply(seq_along(sharedCT), function(X){
  hsaData <- fContent$MTRATIO[fContent$CT %in% sharedCT[X] & fContent$SP %in% 'Homo sapiens']
  mmuData <- fContent$MTRATIO[fContent$CT %in% sharedCT[X] & fContent$SP %in% 'Mus musculus']
  if(length(hsaData) >= 1000 & length(mmuData) >= 1000){
    testOut <- t.test(hsaData,mmuData, conf.level = 0.95)
    c(testOut$estimate, testOut$conf.int, testOut$p.value)  
  } else {
    return(rep(NA, 5))
  }
})
outCT <- t(outCT)
colnames(outCT) <- c('mean Hsa', 'mean Mmu', 'CI low', 'CI high', 'p-value')           
outCT <- as.data.frame(outCT)
outCT$`p-adj` <- p.adjust(outCT$`p-value`, method = 'fdr')
outCT$`Hsa > Mmu` <- outCT$`mean Hsa` > outCT$`mean Mmu`
rownames(outCT) <- sharedCT
outCT <- outCT[complete.cases(outCT),]

      