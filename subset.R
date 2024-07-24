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

project_name ="OCT4subset"
proj_subset <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


#------------------------
#Plotting ATAC Heatmap 
#-------------------------
#-------------------------------
cM_atac_rna <- confusionMatrix(paste0(proj_subset$Clusters_ATAC), paste0(proj_subset$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)

figure_name <- project_name
figure_name <- paste(figure_name, "_atac_rna_heatmap.pdf", sep="")
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
proj_subset <- addImputeWeights(ArchRProj = proj_subset,reducedDims = "LSI_Combined")

saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4subset", load = FALSE)
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
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset),log2Norm=TRUE)
p
dev.off() 



#-----------------------------
#Gex heatmap
#-----------------------------
features <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "GeneExpressionMatrix",
groupBy = "Clusters_Combined",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)
gexHeatmap <- plotMarkerHeatmap(
seMarker = features,
cutOff = "FDR <= 0.05 & Log2FC >= 1",
transpose = TRUE
)
figure_name <- project_name
figure_name <- paste(figure_name,"_gexHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(gexHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
#------------------------
#MarkerGenes heatmap
#------------------------
figure_name = proj_subset
figure_name <- paste(figure_name,"_markersHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
subsetSE <- features[which(rowData(features)$name %in% markerGenes),]
markersHeatmap <- plotMarkerHeatmap(seMarker = subsetSE)
draw(markersHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

featuresControls <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "GeneExpressionMatrix",
groupBy = "Sample",useGroups = "Control_Oct4",
  bgdGroups = "Control_mCherry",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

write.csv(assays(featuresControls)$MeanBGD, "featuresControlsMeanBGD.csv")
write.csv(assays(featuresControls)$MeanDiff, "featuresControlsMeanDiff.csv")
write.csv(assays(featuresControls)$FDR, "featuresControlsFDR.csv")
write.csv(assays(featuresControls)$Log2FC, "featuresControlsLog2FC.csv")
write.csv(rowData(featuresControls)$name,"featuresControlsnames.csv")

featuresRBPJ <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "GeneExpressionMatrix",
groupBy = "Sample",useGroups = "Rbpj_Oct4",
  bgdGroups = "Rbpj_mCherry",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

write.csv(assays(featuresRBPJ)$MeanBGD, "featuresRBPJMeanBGD.csv")
write.csv(assays(featuresRBPJ)$MeanDiff, "featuresRBPJMeanDiff.csv")
write.csv(assays(featuresRBPJ)$FDR, "featuresRBPJFDR.csv")
write.csv(assays(featuresRBPJ)$Log2FC, "featuresRBPJLog2FC.csv")
write.csv(rowData(featuresRBPJ)$name,"featuresRBPJnames.csv")



figure_name = proj_subset
figure_name <- paste(figure_name,"_markersSampleHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
subsetSE <- features[which(rowData(featuresperSample)$name %in% markerGenes),]
markersHeatmap <- plotMarkerHeatmap(seMarker = subsetSE)
draw(markersHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
#--------------------------------
#---------------------------------
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


peaksRBPJ <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "PeakMatrix",
groupBy = "Sample",useGroup = "Rbpj_Oct4",
bgdGroups = "Rbpj_mCherry",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)



peaksControls <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "PeakMatrix",
groupBy = "Sample",useGroups = "Control_Oct4",
bgdGroups = "Control_mCherry",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)



markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_subset,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters_Combined",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_subset,
    useMatrix = "PeakMatrix",
    groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markersPeaks
markersPeaks$C20 

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 2",
  transpose = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#---------------------------
#Plot specific genes' peaks
#---------------------------
p <- plotBrowserTrack(
    ArchRProj = proj_subset,
    groupBy = "Clusters_Combined",
    geneSymbol = c("Pou5f1"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 4", returnGR = TRUE),
    upstream = 10000,
    downstream = 10000
)
grid::grid.newpage()
grid::grid.draw(p$Pou5f1)
figure_name <- project_name
figure_name <- paste(figure_name,"_pou5f1.pdf", sep="")
pdf(file =figure_name, width=12)
dev.off()


p <- plotBrowserTrack(
    ArchRProj = proj_subset,
    groupBy = "Clusters_Combined",
    geneSymbol = c("Gad2"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 4", returnGR = TRUE),
    upstream = 10000,
    downstream = 10000
)
figure_name <- project_name
grid::grid.newpage()
grid::grid.draw(p$Gad2)
figure_name <- paste(figure_name,"_Gad2.pdf", sep="")
pdf(file =figure_name, width=12)
dev.off()



peakSet <- getPeakSet(proj_subset)
write.csv(peakSet, "peakSet.csv")
peakMatrix <- getAvailableMatrices(proj_subset)
peakMatrix


saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4subset", load = FALSE)

#----------------------
#Calling Motifs 
#----------------------
proj_subset <- addMotifAnnotations(ArchRProj = proj_subset, motifSet = "cisbp", name = "Motif")



motifsUPControls <- peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.5 & Log2FC >= 0.5" 
  )

motifsDoControls <-peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.5 & Log2FC >= -0.5",  )


motifsUPRBPJ <- peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.5 & Log2FC >= 0.5"
  )

motifsDoRBPJ <-peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.5 & Log2FC >= -0.5"  )


dfUPControls <-data.frame(TF =rownames(motifsUPControls), mlog10Padj =assay(motifsUPControls)[,1])
dfUPControls <- dfUPControls[order(dfUPControls$mlog10Padj, decreasing = TRUE),]
dfUPControls$rank <-seq_len(nrow(dfUPControls))

write.csv(dfUPControls, "dfUPControls.csv") 

dfUPRBPJ <-data.frame(TF =rownames(motifsUPRBPJ), mlog10Padj =assay(motifsUPRBPJ)[,1])
dfUPRBPJ <- dfUPRBPJ[order(dfUPRBPJ$mlog10Padj, decreasing = TRUE),]
dfUPRBPJ$rank <-seq_len(nrow(dfUPRBPJ))
write.csv(dfUPRBPJ, "dfUPRBPJ.csv") 



heatmapEMUP <- plotEnrichHeatmap(motifsUP, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapEMUP, heatmap_legend_side = "bot", annotation_legend_side = "bot")
figure_name =""
figure_name <- paste(figure_name,"motifsUP.pdf", sep="")
pdf(file =figure_name, width=12)
dev.off()

motifsDo <-peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.5 & Log2FC <= -0.5"  )
heatmapEMDo <- plotEnrichHeatmap(motifsDo, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEMDo, heatmap_legend_side = "bot", annotation_legend_side = "bot")
figure_name =""
figure_name <- paste(figure_name,"motifsDo.pdf", sep="")
pdf(file =figure_name, width=12)
dev.off()

saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4subset", load = FALSE)




#Footprints 


motifs = c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
seFoot <- getFootprints(ArchRProj = proj_subset,positions = motifPositions[markerMotifs],groupBy = "Clusters_Combined")
#figure will be in OCT4subset/Plots 
plotFootprints(seFoot = seFoot,ArchRProj = proj_subset,normMethod = "Subtract",plotName = "Footprints-Subtract-Bias",addDOC = FALSE, smoothWindow = 5)
dev.off() 
