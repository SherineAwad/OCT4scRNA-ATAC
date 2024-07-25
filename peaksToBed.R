library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

peaksfile = paste(filename,".csv",sep="") 
bedfile = paste(filename,".bed", sep="")

df <- readfile <- read.csv(peaksfile, head = TRUE, sep=",")

newdf = df[c("seqnames", "start", "end")]

mybed <- makeGRangesFromDataFrame(newdf,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)

rtracklayer::export.bed(mybed, bedfile, "bed")

