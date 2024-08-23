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

args <- commandArgs(trailingOnly = TRUE)

project_name = args[1]
groupBy ="Sample" 

myProject <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")

featuresControls <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "GeneExpressionMatrix",
groupBy = groupBy, ,useGroups = "Control_Oct4",
  bgdGroups = "Control_mCherry",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

figure_name = project_name
figure_name <- paste(figure_name,"_ControlsmarkersHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
subsetSE <- featuresControls[which(rowData(featuresControls)$name %in% markerGenes),]
ControlsmarkersHeatmap <- plotMarkerHeatmap(seMarker = subsetSE,cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",plotLog2FC = TRUE)
draw(ControlsmarkersHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

featuresRBPJ <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "GeneExpressionMatrix",
groupBy = groupBy,useGroups = "Rbpj_Oct4",
  bgdGroups = "Rbpj_mCherry",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

figure_name = project_name
figure_name <- paste(figure_name,"_RBPJmarkersHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
subsetSE <- featuresRBPJ[which(rowData(featuresRBPJ)$name %in% markerGenes),]
RBPJmarkersHeatmap <- plotMarkerHeatmap(seMarker = subsetSE,cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotLog2FC = TRUE)
draw(RBPJmarkersHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


filename = project_name
filename <- paste(filename, "DEGsControls.csv", sep="_")
df <- data.frame(genes=rowData(featuresControls), Log2FC=assays(featuresControls)$Log2FC, FDR=assays(featuresControls)$FDR, Mean= assays(featuresControls)$Mean, 
MeanDiff=assays(featuresControls)$MeanDiff, MeanBGD=assays(featuresControls)$MeanBGD, Pval=assays(featuresControls)$Pval)
colnames(df) <- c("genes.seqnames", "genes.idx", "genes.start", "genes.end", "genes.name", "genes.strand", "Log2FC","FDR", "Mean", "MeanDiff","MeanBGD", "Pvalue")
write.csv(df, filename,row.names=FALSE)


filename = project_name
filename <- paste(filename, "DEGsRBPJ.csv", sep="_")
df <- data.frame(genes=rowData(featuresRBPJ), Log2FC=assays(featuresRBPJ)$Log2FC, FDR=assays(featuresRBPJ)$FDR, Mean= assays(featuresRBPJ)$Mean,
MeanDiff=assays(featuresRBPJ)$MeanDiff, MeanBGD=assays(featuresRBPJ)$MeanBGD, Pval=assays(featuresRBPJ)$Pval)
colnames(df) <- c("genes.seqnames", "genes.idx", "genes.start", "genes.end", "genes.name", "genes.strand", "Log2FC","FDR", "Mean", "MeanDiff","MeanBGD", "Pvalue")
write.csv(df, filename,row.names=FALSE)

saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)


#================================
#Plot specific genes' peaks
#================================

p <- plotBrowserTrack(
    ArchRProj = myProject,
    groupBy = "Sample",
    geneSymbol = c("Pou5f1"),
    features = getPeakSet(myProject),
    loops = getCoAccessibility(myProject),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Pou5f1)
figure_name <- project_name
figure_name <- paste(figure_name,"_pou5f1.pdf", sep="")
#saved in project Plots folder
plotPDF(p, name = figure_name, width = 5, height = 5, ArchRProj = myProject, addDOC = FALSE)
dev.off()

saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)

