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
proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
#-----------------------------
#Add Clusters
#----------------------------
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.6, force = TRUE)
proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.6, force = TRUE)

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
proj_ALL <- addImputeWeights(ArchRProj = proj_ALL,reducedDims = "LSI_Combined") # scaleDims=TRUE, corCutOff =0.5)
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)
#-------------------------------------
#Plotting Gene Expressions 
#-------------------------------------

markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")


markers1 <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1")
markers2 <- c("Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2") 
markers3 <- c("Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2")
markers4 <- c("Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2") 
markers5 <- c("Csf1r", "Ccr2", "Pax2","Tie1", "Klf4","Grm6") 

figure_name <- project_name
figure_name <- paste(figure_name,"_features.pdf", sep="")
pdf(file =figure_name, width=12)

#Plotting all Genes one gene per page 
p <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markerGenes,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL) ,log2Norm=TRUE)
p
dev.off() 

#Plotting group of genes per page 
p1 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markers1,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL), log2Norm=TRUE)

p2 <-plotEmbedding(
ArchRProj = proj_ALL,
colorBy = "GeneExpressionMatrix",
name = markers2,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_ALL), log2Norm=TRUE)

#Check /nfs/turbo/umms-thahoang/sherine/mouseCutandTag/archr/OCT4RBPJ/Plots for plots 
plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p1 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_ALL,
addDOC = FALSE,
width = 20,
height = 20
)

plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p1 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_ALL,
addDOC = FALSE,
width = 20,
height = 20
)

plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p2 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_ALL,
addDOC = FALSE,
width = 20,
height = 20
)

 


#------------------------------
#Browser Track
#------------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_browserTrack.pdf", sep="")
pdf(file =figure_name, width=12)
> p <- plotBrowserTrack(
ArchRProj = proj_ALL,
groupBy = "Clusters",
geneSymbol = markerGenes,
upstream = 50000,
downstream = 50000
)
dev.off()
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
markerList
markerList$C1 

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


peakSet <- getPeakSet(proj_ALL)
write.csv(peakSet, "peakSet.csv")
peakMatrix <- getAvailableMatrices(proj_ALL)
peakMatrix


saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)

#----------------------
#Calling Motifs 
#----------------------
proj_ALL <- addMotifAnnotations(ArchRProj = proj_ALL, motifSet = "cisbp", name = "Motif")
motifsUP <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_ALL,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
#-----------------------
#Plotting Motifs 
#-----------------------
heatmapEM <- plotEnrichHeatmap(motifsUP, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_ALL, addDOC = FALSE)
dev.off()

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

figure_name =""
figure_name <- paste(figure_name,"motifsUP.pdf", sep="")
pdf(file =figure_name, width=12)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggUp
dev.off()

motifsDo <-peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_ALL,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"  )

motifsDo
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

figure_name =""
figure_name <- paste(figure_name,"motifsDo.pdf", sep="")
pdf(file =figure_name, width=12)

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
dev.off()

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4RBPJ", load = FALSE)



clusters <- c("C3", "C4", "C5", "C6", "C8", "C9", "C10","C12", "C13", "C14", "C15")
proj_subset = subsetArchRProject(proj_ALL,proj_ALL$cellNames[proj_ALL$Clusters_Combined %in% clusters] , outputDirectory = "OCT4_subset" ,force =TRUE,dropCells = TRUE)

}
