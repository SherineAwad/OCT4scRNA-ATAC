#!/bin/bash

#Convert peaks CSV file  to bed file 
Rscript peaksToBed.R OCT4subset_PeaksControls 0
Rscript peaksToBed.R OCT4subset_PeaksRBPJ 0
#Use Hommer annotatePeaks.pl to annotate peaks bed file 
annotatePeaks.pl OCT4subset_PeaksControls.bed mm10 -annStats OCT4subset_PeaksControls.stats > OCT4subset_PeaksControlsAnnotated.bed
annotatePeaks.pl OCT4subset_PeaksRBPJ.bed mm10 -annStats OCT4subset_PeaksRBPJ.stats > OCT4subset_PeaksRBPJAnnotated.bed

#Merge information in peaks CSV file and the annotated peaks file 
python mergePeaksAnnotations.py OCT4subset_PeaksControls.csv OCT4subset_PeaksControlsAnnotated.bed > OCT4subset_PeaksControlsMergedAnnotations.csv
python mergePeaksAnnotations.py OCT4subset_PeaksRBPJ.csv OCT4subset_PeaksRBPJAnnotated.bed > OCT4subset_PeaksRBPJMergedAnnotations.csv


