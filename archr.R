library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)

setwd("/nfs/turbo/umms-thahoang/sherine/mouseCutandTag/archr")
#addArchRThreads(threads = 32) 

addArchRGenome("mm10")


######
atacFiles <- c("Control_mCherry" = "Control_mCherry_NMDA_atac_fragments.tsv.gz", "Control_Oct4" = "Control_Oct4_NMDA_atac_fragments.tsv.gz", "Rbpj_mCherry" = "Rbpj_mCherry_NMDA_atac_fragments.tsv.gz", "Rbpj_Oct4" = "Rbpj_Oct4_NMDA_atac_fragments.tsv.gz")

rnaFiles <- c("Control_mCherry" = "Control_mCherry_NMDA_filtered_feature_bc_matrix.h5", "Control_Oct4" = "Control_Oct4_NMDA_filtered_feature_bc_matrix.h5", "Rbpj_mCherry" = "Rbpj_mCherry_NMDA_filtered_feature_bc_matrix.h5", "Rbpj_Oct4" = "Rbpj_Oct4_NMDA_filtered_feature_bc_matrix.h5")


names(atacFiles)
names(rnaFiles)
all.equal(names(atacFiles), names(rnaFiles))
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

#######

#######
#######
ArrowFiles <- c("Control_mCherry.arrow","Control_Oct4.arrow","Rbpj_mCherry.arrow","Rbpj_Oct4.arrow")
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "OCT4andRBPJ", copyArrows = FALSE)
#need to check indeces of ArrowFiles
proj_A <- ArchRProject(ArrowFiles[1], outputDirectory = "proj_A", copyArrows = FALSE)
proj_B <- ArchRProject(ArrowFiles[2], outputDirectory = "proj_B", copyArrows = FALSE)
proj_C <- ArchRProject(ArrowFiles[3], outputDirectory = "proj_C", copyArrows = FALSE)
proj_D <- ArchRProject(ArrowFiles[4], outputDirectory = "proj_D", copyArrows = FALSE)

seRNA_A <- import10xFeatureMatrix(input="Control_mCherry_NMDA_filtered_feature_bc_matrix.h5", names="Control_mCherry")
seRNA_B <- import10xFeatureMatrix(input="Control_Oct4_NMDA_filtered_feature_bc_matrix.h5", names="Control_Oct4")
seRNA_C <- import10xFeatureMatrix(input="Rbpj_mCherry_NMDA_filtered_feature_bc_matrix.h5", names="Rbpj_mCherry")
seRNA_D <- import10xFeatureMatrix(input="Rbpj_Oct4_NMDA_filtered_feature_bc_matrix.h5", names="Rbpj_Oct4")


gr <- rowRanges(seRNA_A)
seqlevels(gr) <- paste0("chr", seqlevels(gr))
rowRanges(seRNA_A) <- gr
k <- which(seqnames(rowRanges(seRNA_A)) %in% seqnames(proj_A@genomeAnnotation$chromSizes) == T)
seRNA_A = seRNA_A[k,]
k = which(rownames(proj_A@cellColData) %in% colnames(seRNA_A) == T)
proj_A = proj_A[k]


gr <- rowRanges(seRNA_B)
seqlevels(gr) <- paste0("chr", seqlevels(gr))
rowRanges(seRNA_B) <- gr
k <- which(seqnames(rowRanges(seRNA_B)) %in% seqnames(proj_B@genomeAnnotation$chromSizes) == T)
seRNA_B = seRNA_B[k,]
k = which(rownames(proj_B@cellColData) %in% colnames(seRNA_B) == T)
proj_B = proj_B[k]


gr <- rowRanges(seRNA_C)
seqlevels(gr) <- paste0("chr", seqlevels(gr))
rowRanges(seRNA_C) <- gr
k <- which(seqnames(rowRanges(seRNA_C)) %in% seqnames(proj_C@genomeAnnotation$chromSizes) == T)
seRNA_C = seRNA_C[k,]
k = which(rownames(proj_C@cellColData) %in% colnames(seRNA_C) == T)
proj_C = proj_C[k]


gr <- rowRanges(seRNA_D)
seqlevels(gr) <- paste0("chr", seqlevels(gr))
rowRanges(seRNA_D) <- gr
k <- which(seqnames(rowRanges(seRNA_D)) %in% seqnames(proj_D@genomeAnnotation$chromSizes) == T)
seRNA_D = seRNA_D[k,]
k = which(rownames(proj_D@cellColData) %in% colnames(seRNA_D) == T)
proj_D = proj_D[k]

seRNAcombined<-cbind(assay(seRNA_A), assay(seRNA_B),assay(seRNA_C), assay(seRNA_D))
seRNA_all<-SummarizedExperiment(assays=list(counts=seRNAcombined), rowRanges= rowRanges(seRNA_A))

project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "OCT4andRBPJ", copyArrows = FALSE)
proj_ALL <-addGeneExpressionMatrix(input=project_ALL, seRNA=seRNA_all)
proj_ALL <- proj_ALL[proj_ALL$TSSEnrichment > 6 & proj_ALL$nFrags > 2500 & !is.na(proj_ALL$Gex_nUMI)]

table(proj_ALL$Sample)



########
######## 
########
#LSI-ATAC
proj_ALL <- addIterativeLSI(
    ArchRProj = proj_ALL, 
    clusterParams = list(
      resolution = 0.2, 
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix", 
    depthCol = "nFrags",
    name = "LSI_ATAC"
)

##########
#LSI-RNA
proj_ALL <- addIterativeLSI(
    ArchRProj = proj_ALL, 
    clusterParams = list(
      resolution = 0.2, 
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
)

############
proj_ALL <- addCombinedDims(proj_ALL, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)


proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)

proj_ALL <- addGroupCoverages(ArchRProj = proj_ALL, groupBy = "Clusters_Combined")
proj_ALL <- addReproduciblePeakSet(ArchRProj = proj_ALL, groupBy = "Clusters_Combined", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2")
proj_ALL <- addPeakMatrix(ArchRProj = proj_ALL)
proj_ALL <- addPeak2GeneLinks(ArchRProj = proj_ALL, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

p2g <- getPeak2GeneLinks(ArchRProj = proj_ALL)

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)


