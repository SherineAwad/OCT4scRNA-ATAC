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
#addArchRThreads(threads = 32) 

addArchRGenome("mm10")

project_name ="OCT4andRBPJ"

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


if (FALSE)
{
######
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



#-----------------------------------
proj_ALL <- addCombinedDims(proj_ALL, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

#-----------------------------
#Correct for batch effect 
#------------------------------
proj_ALL <-addHarmony(proj_ALL, reducedDims = "LSI_Combined", name = "Harmony", groupBy = "Sample")
proj_ALL <- addUMAP(proj_ALL, reducedDims = "Harmony", name = "UMAP_Harmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")


#-----------------------------
#Add Clusters
#----------------------------
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)
#----------------------------------
#Plotting UMAPS
#----------------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_clusters.pdf", sep="")
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
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)
dev.off()

figure_name <- project_name
figure_name <- paste(figure_name,"_SampleUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2, type = "h")
dev.off()


figure_name <- project_name
figure_name <- paste(figure_name,"_HarmonyUMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")
p2 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Harmony")
ggAlignPlots(p1, p2, type = "h")
dev.off()

#------------------------
#Plotting ATAC Heatmap 
#-------------------------
#-------------------------------
cM_atac_rna <- confusionMatrix(paste0(proj_ALL$Clusters_ATAC), paste0(proj_ALL$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)

figure_name <- project_name
figure_name <- paste(figure_name, "_heatmap.pdf", sep="")
pdf(file =figure_name, width=12)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
dev.off()

#-------------------------------------
#Adding Impute Weights using Harmony
#-------------------------------------
proj_ALL <- addImputeWeights(ArchRProj = proj_ALL,reducedDims = "Harmony", scaleDims=NULL, corCutOff =0.5)

#-------------------------------------
#Plotting Gene Expressions 
#-------------------------------------

markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")
markerGenes1 <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms")
markerGenes2 <- c("Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2")
markerGenes3 <- c("Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6")
markerGenes4 <- c("Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4")
markerGenes5 <- c("Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")



p11 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes1,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_ALL) )

p22 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes2,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_ALL) )


p33 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes3,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_ALL) )

p44 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes4,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_ALL) )


p55 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes5,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_ALL) )

p1 <- lapply(p11, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
figure_name <- project_name
figure_name <- paste(figure_name,"_features1.pdf", sep="")
pdf(file =figure_name, width=12)
do.call(cowplot::plot_grid, c(list(ncol = 3),p1))
dev.off()


p2 <- lapply(p22, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
figure_name <- project_name
figure_name <- paste(figure_name,"_features2.pdf", sep="")
pdf(file =figure_name, width=12)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()


p3 <- lapply(p33, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
figure_name <- project_name
figure_name <- paste(figure_name,"_features3.pdf", sep="")
pdf(file =figure_name, width=12)
do.call(cowplot::plot_grid, c(list(ncol = 3),p3))
dev.off()


p4 <- lapply(p44, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
figure_name <- project_name
figure_name <- paste(figure_name,"_features4.pdf", sep="")
pdf(file =figure_name, width=12)
do.call(cowplot::plot_grid, c(list(ncol = 3),p4))
dev.off()

p5 <- lapply(p55, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
figure_name <- project_name
figure_name <- paste(figure_name,"_features5.pdf", sep="")
pdf(file =figure_name, width=12)
do.call(cowplot::plot_grid, c(list(ncol = 3),p5))
dev.off()

}
#-----------------------------------
#Calling Peaks 
#-----------------------------------

proj_ALL <- addGroupCoverages(ArchRProj = proj_ALL, groupBy = "Clusters_Combined")
proj_ALL <- addReproduciblePeakSet(ArchRProj = proj_ALL, groupBy = "Clusters_Combined", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2")
proj_ALL <- addPeakMatrix(ArchRProj = proj_ALL)
proj_ALL <- addPeak2GeneLinks(ArchRProj = proj_ALL, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = proj_ALL)

#--------------
#Plotting Peaks 
#--------------
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_ALL,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters_Combined",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()

#---------------------------
#Plot specific genes' peaks
#---------------------------
p <- plotBrowserTrack(
    ArchRProj = proj_ALL,
    groupBy = "Clusters_Combined",
    geneSymbol = c("Pou5f1"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Pou5f1)
figure_name <- project_name
figure_name <- paste(figure_name,"_pou5f1.pdf", sep="")
pdf(file =figure_name, width=12)
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)

#----------------------
#Calling Motifs 
#----------------------
proj_ALL <- addMotifAnnotations(ArchRProj = proj_ALL, motifSet = "cisbp", name = "Motif")
motifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_ALL,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
#-----------------------
#Plotting Motifs 
#-----------------------
heatmapEM <- plotEnrichHeatmap(motifs, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)


