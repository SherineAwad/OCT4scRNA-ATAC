library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tidyr)

#This script takes any file with chr in col seqnames, start in start col and end in end col 
#and converts it to bed format 
#It takes 0 or 1 as parameter for 0-base or 1-base genome 
#it takes first parameter file name with no suffix (no .csv for example) and second parameter is 0 or 1
#it outputs same file in parameter plus .bed suffix 

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
base <- as.integer(args[2])

peaksfile = paste(filename,".csv",sep="") 
bedfile = paste(filename,".bed", sep="")

flag = TRUE #default is TRUE

if (base == 1) 
    flag =FALSE
flag 

df <- readfile <- read.csv(peaksfile, head = TRUE, sep=",")

newdf = df[c("seqnames", "start", "end")]

mybed <- makeGRangesFromDataFrame(newdf,ignore.strand=TRUE, starts.in.df.are.0based=flag)

rtracklayer::export.bed(mybed, bedfile, "bed")
