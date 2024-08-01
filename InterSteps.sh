#!/bin/bash

#Convert peaks CSV file  to bed file 
Rscript peaksToBed.R SubsetpeaksControlsCCperSample 0
Rscript peaksToBed.R SubsetpeaksRBPJCCperSample 0 

#Use Hommer annotatePeaks.pl to annotate peaks bed file 
annotatePeaks.pl SubsetpeaksControlsCCperSample.bed mm10 -annStats SubsetpeaksControlsCCperSample.stats > SubsetpeaksControlsCCperSampleAnnotated.bed
annotatePeaks.pl SubsetpeaksRBPJCCperSample.bed mm10 -annStats SubsetpeaksRBPJCCperSample.stats > SubsetpeaksRBPJCCperSampleAnnotated.bed

#Merge information in peaks CSV file and the annotated peaks file 
python mergePeaksAnnotations.py SubsetpeaksControlsCCperSample.csv SubsetpeaksControlsCCperSampleAnnotated.bed > SubsetpeaksControlsCCperSampleMergedAnnotations.csv
python mergePeaksAnnotations.py SubsetpeaksRBPJCCperSample.csv SubsetpeaksRBPJCCperSampleAnnotated.bed > SubsetpeaksRBPJCCperSampleMergedAnnotations.csv


