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

proj_ALL <-addHarmony(
    ArchRProj = proj_ALL,
    reducedDims = "LSI_Combined",
    name = "Harmony",
    groupBy = "Sample")
    
    proj_ALL <- addUMAP(
    ArchRProj = proj_ALL, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine")

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "OCT4andRBPJ", load = FALSE)


