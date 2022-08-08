#File to extract node pairs [edge] for graph construction
#modification for post-grid analysis

import sys

infile = sys.argv[1]
outfile = sys.argv[2]		


ID_seq_dict = {}

masterList = {} #create an empty dictionary

outHandle = open(outfile,'w') 
#ID A-132_20_x6     16      Chromosome10    241402  255     20M     *       0       0       CGTTGCCGAGAGTCGTTTTG    IIIIIIIIIIIIIIIIIIII    XA:i:0  MD:Z:20 NM:i:0

with open(infile) as f:
	for eachLine in f:
		elems = eachLine.split("\t")				# split each element 
		ID = elems[9]		#the sequence is used as the identifier
		chr = elems[2]
		loc = elems[3]
		dir = elems[1]

		strand = "nil"

		if(dir == "0"):
			strand = "+"
		elif(dir == "16"):
			strand = "-"
			
		chr_loc = chr + "_" + loc + "_" + strand
		if(ID in ID_seq_dict):
                        ID_seq_dict[ID].append(chr_loc)
                else:
                        ID_seq_dict[ID] = [chr_loc]


for key, value in ID_seq_dict.iteritems():

	outHandle.write(str(key) + "\t")
	for elems in sorted(value):
		outHandle.write(str(elems) + " , ")
	outHandle.write("\n")

print("processing completed")
	
	
	
