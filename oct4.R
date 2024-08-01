library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(chromVARmotifs)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

setwd("/nfs/turbo/umms-thahoang/sherine/mouseCutandTag/archr")
addArchRThreads(threads = 4) 
addArchRGenome("mm10")
project_name ="OCT4RBPJ"
atacFiles <- c("Control_mCherry" = "Control_mCherry_NMDA_atac_fragments.tsv.gz", "Control_Oct4" = "Control_Oct4_NMDA_atac_fragments.tsv.gz", "Rbpj_mCherry" = "Rbpj_mCherry_NMDA_atac_fragments.tsv.gz", "Rbpj_Oct4" = "Rbpj_Oct4_NMDA_atac_fragments.tsv.gz")

rnaFiles <- c("Control_mCherry" = "Control_mCherry_NMDA_filtered_feature_bc_matrix.h5", "Control_Oct4" = "Control_Oct4_NMDA_filtered_feature_bc_matrix.h5", "Rbpj_mCherry" = "Rbpj_mCherry_NMDA_filtered_feature_bc_matrix.h5", "Rbpj_Oct4" = "Rbpj_Oct4_NMDA_filtered_feature_bc_matrix.h5")


names(atacFiles)
names(rnaFiles)
all.equal(names(atacFiles), names(rnaFiles))
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

#######

#######
#######
ArrowFiles <- c("Control_mCherry.arrow","Control_Oct4.arrow","Rbpj_mCherry.arrow","Rbpj_Oct4.arrow")
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "OCT4RBPJ", copyArrows = FALSE)
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

project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "OCT4RBPJ", copyArrows = FALSE)
proj_ALL <-addGeneExpressionMatrix(input=project_ALL, seRNA=seRNA_all)

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

figure_name <- project_name
figure_name <- paste(figure_name,"_QC.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotGroups(
    ArchRProj = proj_ALL, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
   
p2 <- plotGroups(
    ArchRProj = proj_ALL, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )


p1
p2
dev.off ()

#Before Filtering
table(proj_ALL$Sample)

#Print values we can use for filtering and in plots 
head(proj_ALL@cellColData) 
proj_ALL <- proj_ALL[proj_ALL$TSSEnrichment > 10 & proj_ALL$nFrags > 5000 & !is.na(proj_ALL$Gex_nUMI)]
proj_ALL <- proj_ALL[proj_ALL$Gex_nGenes > 500 & proj_ALL$Gex_nGenes < 5000 & proj_ALL$Gex_nUMI > 1000 & proj_ALL$Gex_nUMI < 15000  ]


#After Filtering 
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

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)

#-----------------------------------
proj_ALL <- addCombinedDims(proj_ALL, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", dimsToUse = 1:20, name = "UMAP_Combined", minDist = 0.4, force = TRUE)

#-----------------------------
#Add Clusters
#----------------------------
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined",dimsToUse = 1:20, name = "Clusters_Combined", resolution = 1.2, force = TRUE)

figure_name <- project_name
figure_name <- paste(figure_name,"_perClustersnUMI.pdf", sep="")
pdf(file =figure_name, width=12)
p <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Clusters_Combined",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p
dev.off()

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)
#----------------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_clustersUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(proj_ALL, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj_ALL, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj_ALL, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
p1 + p2 + p3 + patchwork::plot_layout(nrow = 1, guides = "collect")
p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p))
dev.off()

figure_name <- project_name
figure_name <- paste(figure_name,"_SamplesUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2, type = "h")
dev.off()

#Remove some clusters 

clusters <- c("C2", "C3","C4", "C5", "C6", "C8", "C9", "C10", "C11", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24")
proj_Clean = subsetArchRProject(proj_ALL,proj_ALL$cellNames[proj_ALL$Clusters_Combined %in% clusters] , outputDirectory = "OCT4_Clean" ,force =TRUE,dropCells = TRUE)

proj_Clean <- addUMAP(proj_Clean, reducedDims = "LSI_Combined", dimsToUse = 1:20, name = "UMAP_Combined", minDist = 0.4, force = TRUE)
proj_Clean <- addClusters(proj_Clean, reducedDims = "LSI_Combined",dimsToUse = 1:20, name = "Clusters_Combined", resolution = 1.2, force = TRUE)
proj_Clean <- addImputeWeights(ArchRProj = proj_Clean,reducedDims = "LSI_Combined") #, scaleDims=TRUE, corCutOff =0.5)

figure_name <- project_name
figure_name <- paste(figure_name,"_CleanSamplesUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_Clean, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_Clean, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2, type = "h")
dev.off()

saveArchRProject(ArchRProj = proj_Clean, outputDirectory = "OCT4_Clean", load = FALSE)


Newlabel <- c(
"C1" = "Cone",
"C2" ="Rod",
"C3"= "Cone",
"C4"="Bipolar",
"C5"="Bipolar",
"C6"="Amacrine",
"C7"= "KO MG",
"C8"="KO MG",
"C9"="WT MG",
"C10"="WT MG",
"C11"="WT MG",
"C12"="WT MG",
"C13"="WT MG",
"C14"="KO MG",
"C15"="WT MG",
"C16"="WT MG",
"C17"="KO MG",
"C18"="MGPC",
"C19"="WT MG",
"C20"="KO MG",
"C21"="KO MG",
"C22"="KO MG")


proj_Clean$Celltype <-  mapLabels(proj_Clean$Clusters_Combined, newLabels = Newlabel)


figure_name <- project_name
figure_name <- paste(figure_name,"_NCleanSamplesUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_Clean, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_Clean, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
p3 <- plotEmbedding(ArchRProj = proj_Clean, colorBy = "cellColData", name = "Celltype", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2,p3, type = "h")
dev.off()


saveArchRProject(ArchRProj = proj_Clean, outputDirectory = "OCT4_Clean", load = FALSE)

proj_Clean <- addGroupCoverages(ArchRProj = proj_Clean, groupBy = "Clusters_Combined", force=TRUE)
proj_Clean <- addReproduciblePeakSet(ArchRProj = proj_Clean, groupBy = "Clusters_Combined", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2", force=TRUE)
proj_Clean <- addPeakMatrix(ArchRProj = proj_Clean, force=TRUE)
proj_Clean <- addPeak2GeneLinks(ArchRProj = proj_Clean, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = proj_Clean)
saveArchRProject(ArchRProj = proj_Clean, outputDirectory = "OCT4_Clean", load = FALSE)

#SUBET by celltype WT_MG and KO_MG
celltypes <- c("WT MG","KO MG")
proj_subset = subsetArchRProject(proj_Clean,proj_Clean$cellNames[proj_Clean$Celltype %in% celltypes] , outputDirectory = "OCT4subset" ,force =TRUE,dropCells = TRUE)

