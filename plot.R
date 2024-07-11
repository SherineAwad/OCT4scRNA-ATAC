library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(chromVARmotifs)

setwd("/nfs/turbo/umms-thahoang/sherine/mouseCutandTag/archr")
#addArchRThreads(threads = 32) 

addArchRGenome("mm10")

args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)
#-----------------------
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


#------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_UMAP.pdf", sep="")
pdf(file =figure_name, width=12)
p1 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")
p2 <- plotEmbedding(ArchRProj = proj_ALL, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")
ggAlignPlots(p1, p2, type = "h")
dev.off()
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
#----------------------------
proj_ALL <- addImputeWeights(ArchRProj = proj_ALL,reducedDims = "LSI_Combined", scaleDims=TRUE,corCutOff =0.5)
#-----------------------------
markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")
markerGenes1 <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms") 
markerGenes2 <- c("Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2")
markerGenes3 <- c("Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6")
markerGenes4 <- c("Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4")
markerGenes5 <- c("Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")

#----------------------------
p11 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes1,
quantCut =NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) )
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

p22 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes2,
quantCut =NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) ) 

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


p33 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes3,
quantCut =NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) ) 

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


p44 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes4,
quantCut =NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) ) 

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


p55 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes5,
quantCut =NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) ) 

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
#----------------------
markersGS <-getMarkerFeatures(
    ArchRProj = proj_ALL,
    useMatrix = "GeneExpressionMatrix",
    groupBy = "Clusters_Combined",
    bias =c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
heatmapGS <- plotMarkerHeatmap(                                                                                                                                                                                                                                                                                             seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 1.0",
  labelMarkers = markerGenes,
  transpose = TRUE,
)
figure_name <-project_name
figure_name <- paste(figure_name, "_heatmapfeatures.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()


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

proj_ALL <- addMotifAnnotations(ArchRProj = proj_ALL, motifSet = "cisbp", name = "Motif")


motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_ALL,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(motifsUp, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()


