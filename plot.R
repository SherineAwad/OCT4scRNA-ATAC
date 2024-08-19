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

#Load Clean 

project_name ="MG_Control_subset"
myProject <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)



markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")


#PER CELL TYPE 
######################
#Plot DEGS 
#-----------
DEGs <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "GeneExpressionMatrix",
groupBy = "Celltype",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon"
)

DEGsHeatmap <- plotMarkerHeatmap(
  seMarker = DEGs,
  nPrint = 100,
  clusterCols = TRUE,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  labelMarkers= markerGenes
)


figure_name = project_name
figure_name <- paste(figure_name,"_perCelltypeDEGsHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(DEGsHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot", row_order=c(  "Rod", "Cone", "Amacrine","Bipolar","MGPC", "KO MG","WT MG" ))
dev.off()

#Plot Peaks 
#--------------
peaks <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "PeakMatrix",
groupBy = "Celltype",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)

PeaksHeatmap <- plotMarkerHeatmap(
  seMarker = peaks,
  nPrint = 50,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_perCelltypePeaksHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(PeaksHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot", row_order=c(  "Rod", "Cone", "Amacrine","Bipolar","MGPC", "KO MG","WT MG" ))
dev.off()

#Plot Motifs 
if(FALSE)
{
motifs  <-peakAnnoEnrichment(
    seMarker = peaks,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5")


MotifHeatmap <- plotEnrichHeatmap(motifs, n =30, cutOff=0.5, transpose = TRUE)
figure_name =project_name
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(MotifHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot",row_order=c(  "Rod", "Cone", "Amacrine","Bipolar","MGPC", "KO MG","WT MG" ))
dev.off()
}

saveArchRProject(ArchRProj = myProject, outputDirectory = "OCT4_Clean", load = FALSE)



#Load Subset 
project_name ="OCT4subset"
myProject <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

#PER SAMPLE 
##################

perSampleDEGs <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "GeneExpressionMatrix",
groupBy = "Sample",bias = c("Gex_nUMI","Gex_nGenes"),
testMethod = "wilcoxon",
)

perSampleDEGsHeatmap <- plotMarkerHeatmap(
  seMarker = perSampleDEGs,
  nPrint = 500,
  clusterCols = TRUE,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  pal=paletteContinuous(set = "solarExtra", n = 256, reverse = FALSE),
  transpose = TRUE, plotLog2FC = TRUE,
)

figure_name = project_name
figure_name <- paste(figure_name,"_perSampleDEGsHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(perSampleDEGsHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot") 
dev.off()

#DEGs Results 

Filtered_MG_DEGs <- getMarkers(perSampleDEGs, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
write.csv(Filtered_MG_DEGs, "MGSubsetperSample_DEGs.csv",row.names=FALSE)


perSamplePeaks <- getMarkerFeatures(
ArchRProj = myProject,
useMatrix = "PeakMatrix",
groupBy = "Sample",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)

perSamplePeaksHeatmap <- plotMarkerHeatmap(
  seMarker = perSamplePeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  pal=paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE),
  transpose = TRUE,plotLog2FC = TRUE
)

figure_name <- project_name
figure_name <- paste(figure_name,"_perSamplePeaksHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
draw(perSamplePeaksHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


Filtered_MG_Peaks <- getMarkers(perSamplePeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
write.csv(Filtered_MG_DEGs, "MGSubsetperSample_Peaks.csv",row.names=FALSE)

perSampleMotifs  <-peakAnnoEnrichment(
    seMarker = perSamplePeaks,
    ArchRProj = myProject,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5")




perSampleMotifHeatmap <- plotEnrichHeatmap(perSampleMotifs, n =30, cutOff=0.5, transpose = TRUE)
figure_name =project_name
figure_name <- paste(figure_name,"_perSampleMotifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(perSampleMotifHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


dfMotifs <-data.frame(TF =rownames(perSampleMotifs), mlog10Padj =assay(perSampleMotifs)[,1], mlog10p =assays(perSampleMotifs)[2],  Enrichment =assays(perSampleMotifs)[3], BackgroundProporition= assays(perSampleMotifs)[4], nBackground=assays(perSampleMotifs)[5], BackgroundFrequency=assays(perSampleMotifs)[6], CompareProportion=assays(perSampleMotifs)[7], nCompare=assays(perSampleMotifs)[8],
CompareFrequency=assays(perSampleMotifs)[9],feature=assays(perSampleMotifs)[10] )
write.csv(dfMotifs, "MGSubset_perSampleDFMotifs.csv")



figure_name =project_name
figure_name <- paste(figure_name,"_perSampleEnrichmentMotifs.pdf", sep="")
df = dfMotifs[c("TF", "Enrichment.Control_mCherry", "Enrichment.Control_Oct4")] 
df <- subset(df, Enrichment.Control_mCherry >1.5)
df <- subset(df, Enrichment.Control_Oct4 >1.5)
df$TF <- rownames(df)
df <- df %>% gather(key='sample', value='value', -TF)
ggplot(df, aes(sample, TF)) + geom_tile(aes(fill=value))
dev.off() 

saveArchRProject(ArchRProj = myProject, outputDirectory = "OCT4subset", load = FALSE)

