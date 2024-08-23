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

#==========================================================
#Calling Peaks
#==========================================================

myProject <- addGroupCoverages(ArchRProj = myProject, groupBy = "Clusters_Combined", force=TRUE)
myProject <- addReproduciblePeakSet(ArchRProj = myProject, groupBy = "Clusters_Combined", pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr/bin/macs2", force=TRUE)
myProject <- addPeakMatrix(ArchRProj = myProject, force=TRUE)
myProject <- addPeak2GeneLinks(ArchRProj = myProject, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = myProject)
saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)

#===================================================================
#Plotting Peaks 
#==================================================================
peaksRBPJ <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "PeakMatrix",
groupBy = groupBy,
useGroup = "Rbpj_Oct4",
bgdGroups = "Rbpj_mCherry",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)

peaksControls <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "PeakMatrix",
groupBy = groupBy,
useGroups = "Control_Oct4",
bgdGroups = "Control_mCherry",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)


peaksControlsGR <- getMarkers(
  seMarker = peaksControls,
  cutOff = "FDR >= 0",
  n = NULL,
  returnGR = TRUE
)
filename =project_name 
filename <- paste(filename, "PeaksControls.csv", sep="_")
write.csv(peaksControlsGR, filename,row.names=FALSE)

peaksRBPJGR <- getMarkers(
  seMarker = peaksRBPJ,
  cutOff = "FDR >= 0",
  n = NULL,
  returnGR = TRUE
)
filename =project_name
filename <- paste(filename, "PeaksRBPJ.csv", sep="_")
write.csv(peaksRBPJGR, filename,row.names=FALSE) 

heatmapPeaksControls <- plotMarkerHeatmap(
  seMarker = peaksControls,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksControlsheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaksControls, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


heatmapPeaksRBPJ <- plotMarkerHeatmap(
  seMarker = peaksRBPJ,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_peaksRBPJheatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(heatmapPeaksRBPJ, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)

#==========================================================
#Calling Motifs 
#==========================================================

myProject <- addMotifAnnotations(ArchRProj = myProject, motifSet = "cisbp", name = "Motif", force=TRUE)

motifsUPControls <- peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5" 
  )

motifsDoControls <-peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5")


motifsUPRBPJ <- peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

motifsDoRBPJ <-peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"  )


dfUPControls <-data.frame(TF =rownames(motifsUPControls), mlog10Padj =assay(motifsUPControls)[,1], mlog10p =assays(motifsUPControls)[2],  Enrichment =assays(motifsUPControls)[3], BackgroundProporition= assays(motifsUPControls)[4], nBackground=assays(motifsUPControls)[5], BackgroundFrequency=assays(motifsUPControls)[6], CompareProportion=assays(motifsUPControls)[7], nCompare=assays(motifsUPControls)[8],
CompareFrequency=assays(motifsUPControls)[9],feature=assays(motifsUPControls)[10] ) 
dfUPControls <- dfUPControls[order(dfUPControls$mlog10Padj, decreasing = FALSE),]
dfUPControls$rank <-seq_len(nrow(dfUPControls))
filename =project_name
filename <- paste(filename, "MotifsUPControls.csv", sep="_")
write.csv(dfUPControls, filename) 

dfUPRBPJ <-data.frame(TF =rownames(motifsUPRBPJ), mlog10Padj =assay(motifsUPRBPJ)[,1], mlog10p =assays(motifsUPRBPJ)[2],  Enrichment =assays(motifsUPRBPJ)[3], BackgroundProporition= assays(motifsUPRBPJ)[4], nBackground=assays(motifsUPRBPJ)[5], BackgroundFrequency=assays(motifsUPRBPJ)[6], CompareProportion=assays(motifsUPRBPJ)[7], nCompare=assays(motifsUPRBPJ)[8],
CompareFrequency=assays(motifsUPRBPJ)[9],feature=assays(motifsUPRBPJ)[10] )
dfUPRBPJ <- dfUPRBPJ[order(dfUPRBPJ$mlog10Padj, decreasing = TRUE),]
dfUPRBPJ$rank <-seq_len(nrow(dfUPRBPJ))

filename =project_name
filename <- paste(filename, "MotifsUPRbpj.csv", sep="_")
write.csv(dfUPRBPJ, filename) 


dfDoControls <-data.frame(TF =rownames(motifsDoControls), mlog10Padj =assay(motifsDoControls)[,1], mlog10p =assays(motifsDoControls)[2],  Enrichment =assays(motifsDoControls)[3], BackgroundProporition= assays(motifsDoControls)[4], nBackground=assays(motifsDoControls)[5], BackgroundFrequency=assays(motifsDoControls)[6], CompareProportion=assays(motifsDoControls)[7], nCompare=assays(motifsDoControls)[8],
CompareFrequency=assays(motifsDoControls)[9],feature=assays(motifsDoControls)[10] )
dfDoControls <- dfDoControls[order(dfDoControls$mlog10Padj, decreasing = TRUE),]
dfDoControls$rank <-seq_len(nrow(dfDoControls))
filename =project_name
filename <- paste(filename, "MotifsDoControls.csv", sep="_")
write.csv(dfDoControls, filename)


dfDoRBPJ <-data.frame(TF =rownames(motifsDoRBPJ), mlog10Padj =assay(motifsDoRBPJ)[,1], mlog10p =assays(motifsDoRBPJ)[2],  Enrichment =assays(motifsDoRBPJ)[3], BackgroundProporition= assays(motifsDoRBPJ)[4], nBackground=assays(motifsDoRBPJ)[5], BackgroundFrequency=assays(motifsDoRBPJ)[6], CompareProportion=assays(motifsDoRBPJ)[7], nCompare=assays(motifsDoRBPJ)[8],
CompareFrequency=assays(motifsDoRBPJ)[9],feature=assays(motifsDoRBPJ)[10] )

dfDoRBPJ <- dfDoRBPJ[order(dfDoRBPJ$mlog10Padj, decreasing = TRUE),]
dfDoRBPJ$rank <-seq_len(nrow(dfDoRBPJ))
filename =project_name
filename <- paste(filename, "MotifsDoRbpj.csv", sep="_")
write.csv(dfDoRBPJ, filename)


heatmapUPControls <- plotEnrichHeatmap(motifsUPControls, n = 30, cutOff=0.5, transpose = FALSE)
figure_name = project_name
figure_name <- paste(figure_name,"motifsUPControls.pdf", sep="")
pdf(file =figure_name, width=12) 
ComplexHeatmap::draw(heatmapUPControls, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()



heatmapDoControls <- plotEnrichHeatmap(motifsDoControls, n =30, cutOff=0.5, transpose = FALSE)
figure_name =project_name
figure_name <- paste(figure_name,"motifsDoControls.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapDoControls, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()



heatmapUPRBPJ <- plotEnrichHeatmap(motifsUPRBPJ, n = 30, cutOff=0.5, transpose = FALSE)
figure_name =project_name
figure_name <- paste(figure_name,"motifsUPRBPJ.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapUPRBPJ, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

heatmapDoRBPJ <- plotEnrichHeatmap(motifsDoRBPJ, n = 30, cutOff =0.5, transpose = FALSE)
figure_name =project_name
figure_name <- paste(figure_name,"motifsDoRBPJ.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(heatmapDoRBPJ, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

figure_name =project_name
figure_name <- paste(figure_name,"_MotifControlsUPEnrichmentMotifs.pdf", sep="")

pdf(file =figure_name, width=12)
df = dfUPControls[c("TF", "Enrichment.value")]
df <- subset(df, Enrichment.value >1.5)
df$TF <- rownames(df)
df <- df %>% gather(key='sample', value='value', -TF)
ggplot(df, aes(sample, TF)) + geom_tile(aes(fill=value))
dev.off()

figure_name =project_name
figure_name <- paste(figure_name,"_MotifControlsDoEnrichmentMotifs.pdf", sep="")

pdf(file =figure_name, width=12)
df = dfDoControls[c("TF", "Enrichment.value")]
df <- subset(df, Enrichment.value >1.5)
df$TF <- rownames(df)
df <- df %>% gather(key='sample', value='value', -TF)
ggplot(df, aes(sample, TF)) + geom_tile(aes(fill=value))
dev.off()


figure_name =project_name
figure_name <- paste(figure_name,"_MotifRBPJUPEnrichmentMotifs.pdf", sep="")

pdf(file =figure_name, width=12)
df = dfUPRBPJ[c("TF", "Enrichment.value")]
df <- subset(df, Enrichment.value >1.5)
df$TF <- rownames(df)
df <- df %>% gather(key='sample', value='value', -TF)
ggplot(df, aes(sample, TF)) + geom_tile(aes(fill=value))
dev.off()

figure_name =project_name
figure_name <- paste(figure_name,"_MotifRBPJDoEnrichmentMotifs.pdf", sep="")

pdf(file =figure_name, width=12)
df = dfDoRBPJ[c("TF", "Enrichment.value")]
df <- subset(df, Enrichment.value >1.5)
df$TF <- rownames(df)
df <- df %>% gather(key='sample', value='value', -TF)
ggplot(df, aes(sample, TF)) + geom_tile(aes(fill=value))
dev.off()


saveArchRProject(ArchRProj = myProject, outputDirectory = project_name, load = FALSE)

#======================================
#Plot Footprints for specific Motifs
#======================================

motifs = c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Pou5f1", "Gnat2", "Csf1r", "Sox2",  "Rfx", "Lhx2", "Stat")
motifPositions <- getPositions(myProject)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#seFoot works only if groubBy is by Celltype 
seFoot <- getFootprints(ArchRProj = myProject,positions = motifPositions[markerMotifs],groupBy = "Sample")
#figure will be in OCT4subset/Plots
plot_name =""
plot_name =paste(project_name, "FootprintsSample", sep="_")
plotFootprints(seFoot = seFoot,ArchRProj = myProject,normMethod = "Subtract",plotName = plot_name,addDOC = FALSE, smoothWindow = 5)
