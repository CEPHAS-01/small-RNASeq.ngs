#USAGE
#python filename numFiles KeyFile outFile fileIn_1 fileIn_2 ... fileIn_N

import sys

#read each file from the argument list
numFiles = int(sys.argv[1])
allParams = int(numFiles + 4)

if((len(sys.argv) < 5) or (len(sys.argv) != allParams)):
	print("\nMATRIX FILE:: No Replicate - Error in number of parameters!...exiting script\n\n")
	print"USAGE: python matrixNoRep.py numFiles KeyFile outFile fileIn_1 fileIn_2 ... fileIn_N\n"
	details = '''numFiles = the number of files for which the matrix is to be created
keyFile = the non-redundant list of the keys
outFile = Output filename
fileIn_1 = first file
fileIn_2 = second file
fileIn_N = last file\n'''
	print(details)
else:
	key_file = sys.argv[2]
	out_file = sys.argv[3]
	#open the output file
	outHandle = open(out_file,'w')
	#Open key file and read one line at a time
 	with open(key_file) as kf:
		for eachline in kf:
			temp_list = [0] * numFiles
			kSeq = eachline.strip(' \t\n\r')
			upRange = int(numFiles + 4)

			for i in range(4,upRange):
                		with open(sys.argv[i]) as f:
					for eachline in f:
                        			seq = eachline.strip(' \t\n\r')
						if (kSeq == seq):
                                        		curr = int(temp_list[i-4])
							nw = int(curr + 1)
							temp_list[i-4] = nw
                                		else:
                                			continue

			outHandle.write(str(kSeq) + "\t")
			for ind,item in enumerate(temp_list):
				lastItemIndex = numFiles - 1
				if(ind == lastItemIndex):
					outHandle.write(str(item) + "\n")
				else:
					outHandle.write(str(item) + "\t")
	print("Processing completed!")

