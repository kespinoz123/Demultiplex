#!/usr/bin/env python3.10
import argparse
import numpy as np
import bioinfo
import matplotlib.pyplot as plt

#R2 is index 1 with R1 biological
#R3 is index 2 with R4 biological 

R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" 
R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" 
R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" 
# R1="testR1.fq.gz"
#overwrite for test file

def get_args():
    parser = argparse.ArgumentParser(description="A program for demultiplexing")
    parser.add_argument("-f", "--file", help="File name")
    parser.add_argument("-l", "--length", help="length of file", type=int)
    parser.add_argument("-r", "--records", help="Number of records in file", type=int)
    return parser.parse_args()

args=get_args()

#PS9 on creating a 2D List
#position_array = np.zeros((101,4000000),dtype=int) 
#make it without a define range, so refer back to args
# all_qualscores_array=np.empty([args.length,args.records])
#print(args.length, args.records)

#read as text and not binary
import gzip

with gzip.open(R1, 'rt') as fq:
    index=0
    num_rec=0
    all_qualscores_array=np.empty([args.length,args.records], dtype=np.int8)
    #an integer that takes 16 bytes of storage
    mean = np.empty([args.length])
    for line in fq:
        index+=1
        if index %4 == 0:
            line = line.strip()
            # print(line)
            for position, phred_sc in enumerate(line):  
                qs=bioinfo.convert_phred(phred_sc)
                # print(qs,end="\t")
                all_qualscores_array[position, num_rec]=qs
            num_rec+=1
#this is to check code is working
print(len(all_qualscores_array))
#use len for the length function of X

# print(len(all_qualscores_array[0]))
print(all_qualscores_array)

for index in range(len(all_qualscores_array)):
    mean[index]=np.mean(all_qualscores_array[index])
    print(index, "\t", mean[index])

x =range(args.length)
y= mean

plt.bar(x,y, color="lavender")

plt.ylabel("Mean Quality Score", fontsize=10)
plt.xlabel("Read Position (bp)", fontsize=10)
plt.title("Distribution of Mean Quality Score across Read Position", fontsize=10)
plt.savefig(f"Hist_R1.png")
#plt.savefig(f"Hist_k{args.size}.png")

# Wrapper Notes:
# Its a bash script that calls my python script
# we write a wrapper to pass sbash commands and do conda activate (or anything else we need to do bash wise before we run our python script)
# command to run: sbatch [insert name of script for wrapper]