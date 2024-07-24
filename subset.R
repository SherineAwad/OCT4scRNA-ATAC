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
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markerGenes,
quantCut = NULL,
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset) ,log2Norm=TRUE)
p
dev.off()

#Plotting group of genes per page
p1 <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markers1,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset), log2Norm=TRUE)

p2 <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markers2,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset), log2Norm=TRUE)

p3 <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markers3,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset), log2Norm=TRUE)

p4 <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markers4,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset), log2Norm=TRUE)


p5 <-plotEmbedding(
ArchRProj = proj_subset,
colorBy = "GeneExpressionMatrix",
name = markers5,
quantCut = c(0.01, 0.99),
embedding = "UMAP_Combined",  imputeWeights= getImputeWeights(proj_subset), log2Norm=TRUE)


plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p1 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_subset,
addDOC = FALSE,
width = 20,
height = 20
)

plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p2 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_subset,
addDOC = FALSE,
width = 20,
height = 20
)

plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p3 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_subset,
addDOC = FALSE,
width = 20,
height = 20
)


plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p4 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_subset,
addDOC = FALSE,
width = 20,
height = 20
)

plotPDF(
do.call(cowplot::plot_grid, c(list(ncol = 3), p5 ) ) ,
name = "OCT4RBPJ_features1.pdf",
ArchRProj = proj_subset,
addDOC = FALSE,
width = 20,
height = 20
)


#-----------------------------
#Gex heatmap
#-----------------------------
features <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "GeneExpressionMatrix",
groupBy = "Clusters_Combined", #you can change to Sample
c("Gex_nUMI","Gex_nGenes"),
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


featuresRBPJ <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "GeneExpressionMatrix",
groupBy = "Sample",useGroups = "Rbpj_Oct4",
  bgdGroups = "Rbpj_mCherry",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

df <- data.frame(genes=rowData(featuresControls), Log2FC=assays(featuresControls)$Log2FC, FDR=assays(featuresControls)$FDR, Mean= assays(featuresControls)$Mean, 
MeanDiff=assays(featuresControls)$MeanDiff, MeanBGD=assays(featuresControls)$MeanBGD, Pval=assays(featuresControls)$Pval)


colnames(df) <- c("genes.seqnames", "genes.idx", "genes.start", "genes.end", "genes.name", "genes.strand", "Log2FC","FDR", "Mean", "MeanDiff","MeanBGD", "Pvalue")
write.csv(df, "genes_Controls.csv") 


df <- data.frame(genes=rowData(featuresRBPJ), Log2FC=assays(featuresRBPJ)$Log2FC, FDR=assays(featuresRBPJ)$FDR, Mean= assays(featuresRBPJ)$Mean,
MeanDiff=assays(featuresRBPJ)$MeanDiff, MeanBGD=assays(featuresRBPJ)$MeanBGD, Pval=assays(featuresRBPJ)$Pval)

colnames(df) <- c("genes.seqnames", "genes.idx", "genes.start", "genes.end", "genes.name", "genes.strand", "Log2FC","FDR", "Mean", "MeanDiff","MeanBGD", "Pvalue")


write.csv(df, "genes_RBPJ.csv")

#--------------------------------
#---------------------------------
#-----------------------------------
#Calling Peaks 
#-----------------------------------

proj_subset <- addGroupCoverages(ArchRProj = proj_subset, groupBy = "Sample")
proj_subset <- addReproduciblePeakSet(ArchRProj = proj_subset, groupBy = "Sample", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2")
proj_subset <- addPeakMatrix(ArchRProj = proj_subset)
proj_subset <- addPeak2GeneLinks(ArchRProj = proj_subset, reducedDims = "LSI_ATAC", useMatrix = "GeneExpressionMatrix")
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

heatmapPeaksControls <- plotMarkerHeatmap(
  seMarker = peaksControls,
  cutOff = "FDR <= 0.5 & abs(Log2FC) >= 0.5",
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksControlsheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaksControls, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


heatmapPeaksRBPJ <- plotMarkerHeatmap(
  seMarker = peaksRBPJ,
  cutOff = "FDR <= 0.5 & abs(Log2FC >= 0.5)",
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksRBPJheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaksRBPJ, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#---------------------------
#Plot specific genes' peaks
#---------------------------
p <- plotBrowserTrack(
    ArchRProj = proj_subset,
    groupBy = "Clusters_Combined",
    geneSymbol = c("Pou5f1"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.5 & abs(Log2FC >=0.5", returnGR = TRUE),
    upstream = 10000,
    downstream = 10000
)
grid::grid.newpage()
grid::grid.draw(p$Pou5f1)
figure_name <- project_name
figure_name <- paste(figure_name,"_pou5f1.pdf", sep="")
pdf(file =figure_name, width=12)
dev.off()

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


dfDoControls <-data.frame(TF =rownames(motifsDoControls), mlog10Padj =assay(motifsDoControls)[,1])
dfDoControls <- dfDoControls[order(dfDoControls$mlog10Padj, decreasing = TRUE),]
dfDoControls$rank <-seq_len(nrow(dfDoControls))

write.csv(dfDoControls, "dfDoControls.csv")

dfDoRBPJ <-data.frame(TF =rownames(motifsDoRBPJ), mlog10Padj =assay(motifsDoRBPJ)[,1])
dfDoRBPJ <- dfDoRBPJ[order(dfDoRBPJ$mlog10Padj, decreasing = TRUE),]
dfDoRBPJ$rank <-seq_len(nrow(dfDoRBPJ))
write.csv(dfDoRBPJ, "dfDoRBPJ.csv")


heatmapUPControls <- plotEnrichHeatmap(motifsUPControls, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifsUPControls.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapUPControls, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


heatmapUPRBPJ <- plotEnrichHeatmap(motifsUPRBPJ, n = 7, transpose = TRUE)
figure_name =""
figure_name <- paste(figure_name,"motifsUPRBPJ.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapUPRBPJ, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveArchRProject(ArchRProj = proj_subset, outputDirectory = "OCT4subset", load = FALSE)


motifs = c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
seFoot <- getFootprints(ArchRProj = proj_subset,positions = motifPositions[markerMotifs],groupBy = "Clusters_Combined")
#figure will be in OCT4subset/Plots 
plotFootprints(seFoot = seFoot,ArchRProj = proj_subset,normMethod = "Subtract",plotName = "Footprints-Subtract-Bias",addDOC = FALSE, smoothWindow = 5)
dev.off() 
