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

project_name ="OCT4_subset"
proj_subset <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


#------------------------
#Plotting ATAC Heatmap 
#-------------------------
#-------------------------------
cM_atac_rna <- confusionMatrix(paste0(proj_subset$Clusters_ATAC), paste0(proj_subset$Clusters_RNA))
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
proj_subset <- addImputeWeights(ArchRProj = proj_subset,reducedDims = "Harmony", scaleDims=TRUE)

saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4_subset", load = FALSE)
#-------------------------------------
#Plotting Gene Expressions 
#-------------------------------------
figure_name <- project_name
figure_name <- paste(figure_name,"_features.pdf", sep="")
markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")


p <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markerGenes,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Harmony",  imputeWeights= getImputeWeights(proj_subset) )

p
dev.off() 


#-----------------------------------
#Calling Peaks 
#-----------------------------------

proj_subset <- addGroupCoverages(ArchRProj = proj_subset, groupBy = "Clusters_Combined")
proj_subset <- addReproduciblePeakSet(ArchRProj = proj_subset, groupBy = "Clusters_Combined", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2")
proj_subset <- addPeakMatrix(ArchRProj = proj_subset)
proj_subset <- addPeak2GeneLinks(ArchRProj = proj_subset, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = proj_subset)

#--------------
#Plotting Peaks 
#--------------
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_subset,
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
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_subset, addDOC = FALSE)
dev.off()

#---------------------------
#Plot specific genes' peaks
#---------------------------
p <- plotBrowserTrack(
    ArchRProj = proj_subset,
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
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj_subset, addDOC = FALSE)
dev.off()


peakSet <- getPeakSet(proj_subset)
write.csv(peakSet, "peakSet.csv")
peakMatrix <- getAvailableMatrices(proj_subset)
peakMatrix


saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4andRBPJ", load = FALSE)

#----------------------
#Calling Motifs 
#----------------------
proj_subset <- addMotifAnnotations(ArchRProj = proj_subset, motifSet = "cisbp", name = "Motif")
motifsUP <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_subset,
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
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_subset, addDOC = FALSE)
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
    ArchRProj = proj_subset,
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

saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4andRBPJ", load = FALSE)



clusters <- c("C3", "C4", "C5", "C6", "C8", "C9", "C10","C12", "C13", "C14", "C15")
proj_subset = subsetArchRProject(proj_subset,proj_subset$cellNames[proj_subset$Clusters_Combined %in% clusters] , outputDirectory = "OCT4_subset" ,force =TRUE,dropCells = TRUE)


