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
groupBy =args[2] 

proj_subset <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

#Plotting Gene Expressions 
#-------------------------------------
markerGenes <- c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Sebox","Pou5f1", "Gnat2", "Csf1r")


featuresControls <- getMarkerFeatures(
ArchRProj = proj_subset,
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
ArchRProj = proj_subset,
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

saveArchRProject(ArchRProj = proj_subset, outputDirectory = project_name, load = FALSE)

#--------------------------------
#--------------------------------
#--------------------------------
#Plotting Peaks 
#--------------
peaksRBPJ <- getMarkerFeatures(
ArchRProj = proj_subset,
useMatrix = "PeakMatrix",
groupBy = groupBy,
useGroup = "Rbpj_Oct4",
bgdGroups = "Rbpj_mCherry",
bias = c("TSSEnrichment", "log10(nFrags)"),
testMethod = "wilcoxon"
)

peaksControls <- getMarkerFeatures(
ArchRProj = proj_subset,
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

#----------------------------
#---------------------------
#Plot specific genes' peaks
#---------------------------
p <- plotBrowserTrack(
    ArchRProj = proj_subset,
    groupBy = "Sample",
    geneSymbol = c("Pou5f1"),
    features = getPeakSet(proj_subset),
    loops = getCoAccessibility(proj_subset),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Pou5f1)
figure_name <- project_name
figure_name <- paste(figure_name,"_pou5f1.pdf", sep="")
#saved in project Plots folder
plotPDF(p, name = figure_name, width = 5, height = 5, ArchRProj = proj_subset, addDOC = FALSE)
dev.off()

saveArchRProject(ArchRProj = proj_subset, outputDirectory = project_name, load = FALSE)

#----------------------
#Calling Motifs 
#----------------------
proj_subset <- addMotifAnnotations(ArchRProj = proj_subset, motifSet = "cisbp", name = "Motif", force=TRUE)

motifs  <-peakAnnoEnrichment(
    seMarker = peaks,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR >0")

motifsHeatmap <- plotEnrichHeatmap(motifs, n = 20, cutOff=0.5, binaryClusterRows=FALSE)
figure_name = project_name
figure_name <- paste(figure_name,"motifs.pdf", sep="")
pdf(file =figure_name, width=12)
ComplexHeatmap::draw(motifsHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


motifsUPControls <- peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5" 
  )

motifsDoControls <-peakAnnoEnrichment(
    seMarker = peaksControls,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5")


motifsUPRBPJ <- peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = proj_subset,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

motifsDoRBPJ <-peakAnnoEnrichment(
    seMarker = peaksRBPJ,
    ArchRProj = proj_subset,
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


saveArchRProject(ArchRProj = proj_subset, outputDirectory = project_name, load = FALSE)


motifs = c("Rbfox3", "Sebox", "Gad1", "Elavl3","Sox9", "Glul","Pou4f2", "Rbpms", "Lhx1","Csf1r", "Ccr2", "Pax2","Kcnj8","Rlbp1", "Ascl1", "Otx2", "Olig2", "Crx","Neurog2","Rpe65", "Acta2", "Tie1", "Klf4","Grm6","Grik1","Rho", "Arr3", "Tfap2b", "Vsx1","Insm1","Prdm1", "Elavl4","Gnat1", "Pcp2", "Prkca","Cabp5","Isl1","Slc6a9","Gad2","Chat","Pou5f1", "Gnat2", "Csf1r", "Sox2",  "Rfx", "Lhx2", "Stat")
motifPositions <- getPositions(proj_subset)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#seFoot works only if groubBy is by Celltype 
seFoot <- getFootprints(ArchRProj = proj_subset,positions = motifPositions[markerMotifs],groupBy = "Celltype")
#figure will be in OCT4subset/Plots
plot_name =""
plot_name =paste(project_name, "FootprintsCelltype", sep="_")
plotFootprints(seFoot = seFoot,ArchRProj = proj_subset,normMethod = "Subtract",plotName = plot_name,addDOC = FALSE, smoothWindow = 5)
