library(XML)
library(xml2)
library(pbapply)
library(Matrix)
library(curl)

mainPage <- "https://panglaodb.se/samples.html?species=both&protocol=all%20protocols"
download_html(mainPage, file = "mainPage.txt")
mainPage <- readHTMLTable("mainPage.txt")[[1]]
mainPage <- mainPage[mainPage$`Tissue/Site ` != "",]
mainPage <- as.data.frame.array(mainPage)
mainPage$`Tissue/Site ` <- gsub('\\?','',mainPage$`Tissue/Site `)

dir.create('Results/')
MT <- lapply(seq_len(nrow(mainPage)), function(Z){
  message(Z)
  if(!file.exists(paste0("scDB/",mainPage[Z,3],"_",mainPage[Z,4],".csv"))){
    download.file(paste0("https://panglaodb.se/data_dl.php?sra=",mainPage[Z,2],"&srs=",mainPage[Z,3],"&filetype=R&datatype=readcounts"), destfile = "countData.rdata", method = "curl")
    try(load("countData.rdata"), silent = TRUE)
    if('sm' %in% ls()){
      
      clustersPage <- paste0("https://panglaodb.se/data_dl.php?sra=",mainPage[Z,2],"&srs=",mainPage[Z,3],"&datatype=clusters&filetype=txt")
      download_html(clustersPage, "clusterID.txt")
      clusterID <- read.table("clusterID.txt", header = FALSE, stringsAsFactors = FALSE, row.names = 1)
      
      clusterNames <- paste0("https://panglaodb.se/list_clusters_and_cell_types.html?sra=",mainPage[Z,2],"&srs=",mainPage[Z,3])
      download_html(clusterNames, "clusterNames.txt")
      clusterNames <- readHTMLTable("clusterNames.txt")[[1]]
      clusterNames <- as.data.frame.array(clusterNames)
      
      clusterID <- clusterID[clusterID$V2 %in% clusterNames$`Cluster ID`,, drop = FALSE]
      cID <- as.factor(clusterID$V2)
      levels(cID) <- as.factor(clusterNames$`Inferred cell type`)
      
      sm <- sm[,rownames(clusterID)]
      
      cellType <- cID
      SRA <- mainPage$`SRA `[Z]
      SRS <- mainPage$`SRS `[Z]
      Protocol <- mainPage$`Protocol `[Z]
      SP <- mainPage$`Species `[Z]
      Tissue <- mainPage$`Tissue/Site `[Z]
      cellNames <- colnames(sm)
      librarySize <- colSums(sm)
      nGenes <- colSums(sm != 0)
      geneNames <- rownames(sm)
      mtCounts <- colSums(sm[grepl('^MT-',geneNames, ignore.case = TRUE),])
      mtRatio <- mtCounts/librarySize
      
      
      mtC <- data.frame(SRA = SRA, SRS = SRS, SP = SP, PROTOCOL = Protocol, TISSUE = Tissue, CELL = cellNames, NGENES = nGenes, TOTALCOUNT = librarySize, MTCOUNT = mtCounts, MTRATIO = mtRatio, CT = cellType)
      write.csv(mtC, paste0("Results/",mainPage[Z,3],"_",mainPage[Z,4],".csv"))
    }
  }
})
