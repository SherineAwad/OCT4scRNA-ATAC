#! /usr/bin/env python
import sys
import argparse
import math


def readAnnotations(annofile): 
    peaks ={} 
    for line in open(annofile):
        records = (line.strip()).split("\t")
        peakID = records[1]+":"+records[2]+"-"+records[3]
        peakRecord =""
        for i in range(4,len(records)) :
              peakRecord+=records[i]+"\t" 
        peaks[peakID] =peakRecord
    return (peaks)

def readPeakCSV(csvfile, peaks):
    print("Chr:start-end\twidth\tLog2FC\tFDR\tMeanDif\tStrand\tPeakScore\tFocusRatio/RegionSize\tAnnotation\tDetailed Annotation\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl \tGene Name\tGene Alias\tGene Description\tGene Type")
    for line in open(csvfile):
        if "group_name" in line: 
            continue #Skip header
        records = (line.strip()).split(",")
        start = int(records[3])+ 1
        end = records[4]
        ch = records[2].strip('"')
        peakID =ch+":"+str(start)+"-"+end
        peakRecord =""
        for i in range(5,len(records)) :
              if (i==6):
                  continue
              peakRecord+=records[i]+"\t"
        if peakID in peaks:
            print(peakID+"\t"+peakRecord+peaks[peakID])
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csvfile')
    parser.add_argument('annofile') 
    args = parser.parse_args()
    annofile = args.annofile  
    csvfile = args.csvfile
    peaks = readAnnotations(annofile) 
    readPeakCSV(csvfile, peaks)    
if __name__ == '__main__':
    main()

   

