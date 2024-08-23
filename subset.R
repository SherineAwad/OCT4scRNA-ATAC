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
proj_subset <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)


clusters <- c("C2", "C3","C4", "C5", "C6", "C8", "C9", "C10", "C11", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24")
proj_Clean = subsetArchRProject(proj_ALL,proj_ALL$cellNames[proj_ALL$Clusters_Combined %in% clusters] , outputDirectory = "OCT4_Clean" ,force =TRUE,dropCells = TRUE)
saveArchRProject(ArchRProj = proj_Clean, outputDirectory = "OCT4_Clean", load = FALSE)




Newlabel <- c(
"C1" = "Cone",
"C2" ="Rod",
"C3"= "Cone",
"C4"="Bipolar",
"C5"="Bipolar",
"C6"="Amacrine",
"C7"= "KO MG",
"C8"="KO MG",
"C9"="WT MG",
"C10"="WT MG",
"C11"="WT MG",
"C12"="WT MG",
"C13"="WT MG",
"C14"="KO MG",
"C15"="WT MG",
"C16"="WT MG",
"C17"="KO MG",
"C18"="MGPC",
"C19"="WT MG",
"C20"="KO MG",
"C21"="KO MG",
"C22"="KO MG")


proj_Clean$Celltype <-  mapLabels(proj_Clean$Clusters_Combined, newLabels = Newlabel)

saveArchRProject(ArchRProj = proj_Clean, outputDirectory = "OCT4_Clean", load = FALSE)


#SUBET by celltype WT_MG and KO_MG
celltypes <- c("WT MG","KO MG")
proj_subset = subsetArchRProject(proj_Clean,proj_Clean$cellNames[proj_Clean$Celltype %in% celltypes] , outputDirectory = "OCT4subset" ,force =TRUE,dropCells = TRUE)

