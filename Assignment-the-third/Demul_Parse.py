#!/bin/env python3.10

import gzip
import argparse
import itertools
import bioinfo

def get_args():
	parser = argparse.ArgumentParser(description="A program to parse demultiplex files")
	parser.add_argument("-R1", help="File name", required=True)
	parser.add_argument("-R2", help="File name",required=True)
	parser.add_argument("-R3", help="File name", required=True)
	parser.add_argument("-R4", help="File name", required=True)
	parser.add_argument("-i", help="Index Set", required=True)

	return parser.parse_args()
args=get_args()

def rev_comp(DNA):
    comp = ""
    for base in DNA:
        if base=="A":
            comp+="T"
        elif base=="T":
            comp+="A"
        elif base=="G":
            comp+="C"
        elif base=="C":
            comp+="G"
        elif base=="N":
            comp+="N"
        else:
            print('Input is not a DNA base')
    return comp[::-1]

#--------orignal 4 files to use at the end ------------------
# R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
# R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" 
# R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" 
# R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" 

#---------- original index files to use at the end ------------
# index_file="/projects/bgmp/shared/2017_sequencing/indexes.txt"

#---------- Test Files to use first ------------
# R1="Test_Input_R1.fq.gz"
# R2="Test_Input_R2.fq.gz"
# R3="Test_Input_R3.fq.gz"
# R4="Test_Input_R4.fq.gz"

# index_file="test_indexes.txt"

#---------- Use ARGS ------------
R1=args.R1
R2=args.R2
R3=args.R3
R4=args.R4

index_file=args.i

#---------- Dictionaries to count each condition ------------
matched_counts={"Matched:": 0}
hopped_counts={"Hopped:": 0}
unknown_counts={"Unknown:": 0}

#---------- Opening Index File ------------
index_set=set()

with open(index_file, 'rt') as index_file_handle:
	for line in index_file_handle:
		index=line.strip().split()[4]
	
		index_set.add(index) 
	# print(index_set, sep="\n")
print(index_set)


#-------- Opening all 52 files -------------------
barcodes={}
for barcode in index_set: 
	matched_R1= open(barcode+'_R1.fq', 'w')
	matched_R2= open(barcode+'_R2.fq', 'w')
	barcodes[barcode]=[matched_R1, matched_R2]

#---------- Opening Hopped and Unknown Files for R1 and R2------------
Hopped_list = [open('Hopped_R1.fq', 'w'), open('Hopped_R2.fq', 'w')]

Unknown_list = [open('Unknown_R1.fq', 'w'), open('Unknown_R2.fq', 'w')]
	

#---------- Using Itertools to count each pair product ------------
all_pairs={}
products=itertools.product(index_set, index_set)

print(products)
for pair in products:
	all_pairs[pair]=0
	# print(all_pairs)


#---------- Opening each file and Reading through each Line: Header, Sequence, + and QS ------------
with gzip.open(R1, 'rt') as fq1, gzip.open(R2, 'rt') as fq2, gzip.open(R3, 'rt') as fq3, gzip.open(R4, 'rt') as fq4:
	while True:
		
		header_fq1= fq1.readline().strip()
		if header_fq1=="":
			break
		sequence_fq1= fq1.readline().strip()
		plus_sign_fq1= fq1.readline().strip()
		qs_fq1= fq1.readline().strip()

		header_fq2= fq2.readline().strip()
		sequence_fq2= fq2.readline().strip()
		plus_sign_fq2= fq2.readline().strip()
		qs_fq2= fq2.readline().strip()

		header_fq3= fq3.readline().strip()
		sequence_fq3= fq3.readline().strip()
		plus_sign_fq3= fq3.readline().strip()
		qs_fq3= fq3.readline().strip()

		header_fq4= fq4.readline().strip()
		sequence_fq4= fq4.readline().strip()
		plus_sign_fq4= fq4.readline().strip()
		qs_fq4= fq4.readline().strip()
		
		# print(header_fq1, header_fq2, header_fq3, header_fq4, sep="\n")
		# print("----------")
		
		written = False
#---------- IF statements for each condition: Unknown, Matched and Hopped ------------
		if sequence_fq2 not in index_set or rev_comp(sequence_fq3) not in index_set:
			unknown_counts["Unknown:"]+= 1
			#print('Unknown Match for R2:', sequence_fq2,"\n",'Unknown Match for R3:',rev_comp(sequence_fq3))
			#---------- Writing out Unknown Files  ------------
			Unknown_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
			Unknown_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
		
		elif sequence_fq2 == rev_comp(sequence_fq3):
			for phred_sc in qs_fq2:
				if bioinfo.convert_phred(phred_sc) <30:
					#print('Unknown Match - [ QS < 30 ]:', sequence_fq2, rev_comp(sequence_fq3))
					unknown_counts["Unknown:"]+= 1
					Unknown_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
					Unknown_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
					written = True
					break
			#---------- If statement to avoid double counting in the FWD direction and RVS direction  ------------
			if written == False:
				for phred_sc in qs_fq3:
					if bioinfo.convert_phred(phred_sc) <30:
						#print('Unknown Match - [ QS < 30 ]:', sequence_fq2, rev_comp(sequence_fq3))
						unknown_counts["Unknown:"]+= 1
						Unknown_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
						Unknown_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
						written = True
						break
				#---------- Writing out Dual Matched Files  ------------
				if written == False:
					
					matched_counts["Matched:"]+=1
					#print('Dual Match for R2:', sequence_fq2,"\n",'Dual Match for R3:', rev_comp(sequence_fq3))
					barcodes[sequence_fq2][0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
					barcodes[sequence_fq2][1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
					
					all_pairs[(sequence_fq2,rev_comp(sequence_fq3))]+=1
			
		else:
			for phred_sc in qs_fq2:
				if bioinfo.convert_phred(phred_sc) <30:
					#print('Unknown Match - [ QS < 30 ]:', sequence_fq2, rev_comp(sequence_fq3))
					unknown_counts["Unknown:"]+= 1
					Unknown_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
					Unknown_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
					written = True
					break
			if written == False:
				for phred_sc in qs_fq3:
					if bioinfo.convert_phred(phred_sc) <30:
						#print('Unknown Match - [ QS < 30 ]:', sequence_fq2, rev_comp(sequence_fq3))
						unknown_counts["Unknown:"]+= 1
						Unknown_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
						Unknown_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
						written = True
						break
			#---------- Writing out Dual Hopped Files  ------------
			if written == False:
				hopped_counts["Hopped:"]+=1
				#print('Hopped Match for R2: ', sequence_fq2,"\n", 'Hopped Match for R3: ', rev_comp(sequence_fq3))
				Hopped_list[0].write(header_fq1+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq1+"\n"+plus_sign_fq1+"\n"+qs_fq1+"\n")
				Hopped_list[1].write(header_fq4+"-"+sequence_fq2+"-"+rev_comp(sequence_fq3)+"\n"+sequence_fq4+"\n"+plus_sign_fq4+"\n"+qs_fq4+"\n")
				#print('Hopped Match R2: ', sequence_fq2, 'Hopped Match R3: 'rev_comp(sequence_fq3))
				all_pairs[(sequence_fq2,rev_comp(sequence_fq3))]+=1
				written = True

# -------------- Verifying that my counts for each conditions works -------------- 
#Just print out the number for each condition being tab separated 
print('Total Dual Match: ', matched_counts, 'Total Hopped: ', hopped_counts,'Total Unknown: ', unknown_counts, sep="\t")


#------------ Closing all 52 files------------ 

Unknown_list[0].close()
Unknown_list[1].close()

Hopped_list[0].close()
Hopped_list[1].close()


for barcode in barcodes:
	barcodes[barcode][0].close()
	barcodes[barcode][1].close()

#------------ STATS ----------
#Calculating the sum for each condition and overall sum of all indexes

Total_Matched = sum(matched_counts.values())
print("Total of Dual Matched Indexes", Total_Matched)

Total_Hopped = sum(hopped_counts.values())
print("Total of Hopped Indexes", Total_Hopped)

Total_Unknown = sum(unknown_counts.values())
print("Total of Unknown Indexes", Total_Unknown)

Total_Sum = Total_Matched + Total_Hopped + Total_Unknown
print("Sum of all Index Hits", Total_Sum)


with open("Statistics.md", "wt") as stats:
	stats.write("----------		SUMMARY TABLE	----------"+"\n")
	stats.write("\n")

	stats.write("Dual Matches: " + str(Total_Matched)+"\n")
	stats.write("Hopped Indexes: " + str(Total_Hopped)+"\n")
	stats.write("Unknown Indexes: " + str(Total_Unknown)+"\n")
	stats.write("Total Index Hits in all Files: " + str(Total_Sum)+"\n")

	stats.write("\n")
	stats.write("----------		PERCENTAGES TABLE	----------"+"\n")
	stats.write("Dual Match % = " + str((Total_Matched/Total_Sum)*100)+"%"+"\n")	
	stats.write("Hopped % = " + str((Total_Hopped/Total_Sum)*100)+"%"+"\n")
	stats.write("Unknown % = " + str((Total_Unknown/Total_Sum)*100)+"%"+"\n")
	stats.write("\n")
	stats.write("----------		ALL POSSIBLE PAIR PRODUCTS		----------"+"\n")
	
	for pairs in all_pairs:
		stats.write(pairs[0]+"-"+pairs[1]+": "+str(all_pairs[pairs])+"\n")

# ---------- Calculating the Percentages of each index ------------
# knowing the total will help calculate the percentage of each ex: (Match/Total)*100 = Percentage of Dual Match Indexes
Matched_Percent=(Total_Matched/Total_Sum)*100
print("Dual Matched % ", Matched_Percent)

Hopped_Percent=(Total_Hopped/Total_Sum)*100
print("Index Hopped % ", Hopped_Percent)

Unknown_Percent=(Total_Unknown/Total_Sum)*100
print("Unknown Index % ", Unknown_Percent)


