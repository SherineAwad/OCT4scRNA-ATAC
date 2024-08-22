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
proj_subset <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

#Add ImputeWeights if not added already 
proj_subset <- addImputeWeights(proj_subset,reducedDims="LSI_Combined")

#Trajectory 
trajectory <- c("WT MG", "KO MG", "MGPC", "Bipolar") 
proj_subset <- addTrajectory(
    ArchRProj = proj_subset, 
    name = "trajectory", 
    groupBy = "Celltype",
    trajectory = trajectory, 
    embedding = "UMAP_Combined", 
    force = TRUE
)
p <- plotTrajectory(proj_subset, trajectory = "trajectory", colorBy = "cellColData", embedding="UMAP_Combined", imputeWeights = 
    getImputeWeights(proj_subset),name ="trajectory", continuousSet = "blueYellow") 

p[[1]]

plotPDF(p, name = "TrajectoryUMAP.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 5, height = 5)


trajGIM <- getTrajectory(ArchRProj = proj_subset, name = "trajectory", useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)



#ADD Motif Matrix - Motif deviation 
if("Motif" %ni% names(proj_subset@peakAnnotation)){
    proj_subset <- addMotifAnnotations(ArchRProj = proj_subset, motifSet = "cisbp", name = "Motif")
}


proj_subset <- addBgdPeaks(proj_subset)

proj_subset <- addDeviationsMatrix(
  ArchRProj = proj_subset, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj_subset, name = "MotifMatrix", plot = TRUE)

plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_subset, addDOC = FALSE)

saveArchRProject(ArchRProj = proj_Clean, outputDirectory = project_name, load = FALSE)

