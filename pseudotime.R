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
myProject <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

#Add ImputeWeights if not added already 
myProject <- addImputeWeights(myProject,reducedDims="LSI_Combined")

#Trajectory 
trajectory <- c("WT MG", "KO MG", "MGPC", "Bipolar") 
myProject <- addTrajectory(
    ArchRProj = myProject, 
    name = "trajectory", 
    groupBy = "Celltype",
    trajectory = trajectory, 
    embedding = "UMAP_Combined", 
    force = TRUE
)
p <- plotTrajectory(myProject, trajectory = "trajectory", colorBy = "cellColData", embedding="UMAP_Combined", imputeWeights = 
    getImputeWeights(myProject),name ="trajectory", continuousSet = "blueYellow") 
p
plotPDF(p, name = "TrajectoryUMAP.pdf", ArchRProj = myProject, addDOC = FALSE, width = 5, height = 5)

#Change Matrix to peaks or motif matrix, however motif matrix will need to run the next steps first 
trajGIM <- getTrajectory(ArchRProj = myProject, name = "trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)
trajPeaks <- getTrajectory(ArchRProj = myProject, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = FALSE)

#You can change returnMatrix parameter to True if you need to print the matrix of the heatmap 
p1 <- plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "solarExtra"))
p2 <- plotTrajectoryHeatmap(trajPeaks, pal = paletteContinuous(set = "solarExtra"))

plotPDF(p1, p2,  name = "Traj-Heatmaps.pdf", ArchRProj = myProject, addDOC = FALSE, width = 6, height = 8)


#ADD Motif Matrix - Motif deviation 
if("Motif" %ni% names(myProject@peakAnnotation)){
    myProject <- addMotifAnnotations(ArchRProj = myProject, motifSet = "cisbp", name = "Motif")}



myProject <- addBgdPeaks(myProject)

myProject <- addDeviationsMatrix(
  ArchRProj = myProject, 
  peakAnnotation = "Motif",
  force = FALSE
)

plotVarDev <- getVarDeviations(myProject, name = "MotifMatrix", plot = TRUE)

plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = myProject, addDOC = FALSE)


motifs <- c("Pou5f1", "Rfx4", "Klf4")
markerMotifs <- getFeatures(myProject, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

figure_name =project_name
figure_name <- paste(figure_name,"_MotifUMAP.pdf", sep="")
pdf(file =figure_name, width=12)

#You can use imputeWeights=getImputeWeights etc 
p <- plotEmbedding(
    ArchRProj = myProject,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding ="UMAP_Combined",
    imputeWeights = NULL)
p
dev.off()

saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)

