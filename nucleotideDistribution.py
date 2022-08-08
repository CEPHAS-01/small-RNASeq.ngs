import sys

inFile = sys.argv[1]
outFile = sys.argv[2]
op = int(sys.argv[3])
minLen = int(sys.argv[4])
maxLen = int(sys.argv[5])
# python nucleotideDistribution.py infile outfile operation(1-dict, 2-list) minLen maxLen

outHandle = open(outFile, 'w')
outFileP = "percentage-" + outFile
outHandleP = open(outFileP, 'w')

outFilePlot = "percentage-Plot-" + outFile
outHandlePlot = open(outFilePlot, 'w')


def placeNucleotide(tempList, nucleotide):
	if(nucleotide[:1] == "A"):
		curValue = tempList[0]
		newValue = curValue + 1
		tempList[0] = newValue
	elif(nucleotide[:1] == "T"):
		curValue = tempList[1]
                newValue = curValue + 1
                tempList[1] = newValue
	elif(nucleotide[:1] == "C"):
		curValue = tempList[2]
                newValue = curValue + 1
                tempList[2] = newValue
	elif(nucleotide[:1] == "G"):
		curValue = tempList[3]
                newValue = curValue + 1
                tempList[3] = newValue

	return tempList

def getNucleotideDistro(file, op, minLength, maxlength):
	outDict = {}
	with open(file) as inF:
		for line in inF:
			elem = line.strip(' \t\r\n')
			if((len(elem) >= minLength) and (len(elem) <= maxlength)):
				seqLen = len(elem)
				if(seqLen in outDict):
					oldVal = outDict[seqLen]
					newVal = oldVal + 1
					outDict[seqLen] = newVal
				elif(seqLen not in outDict):
					outDict[seqLen] = 1

	for k,v in sorted(outDict.iteritems()):
		print(str(k) + "\t" + str(v) + "\n")


def getNucleotidePosition(file, op, minLength, maxlength):
	outDict = {}
	outList = [0,0,0,0]
        with open(file) as inF:
                for line in inF:
			elem = line.strip(' \t\r\n')
			if((len(elem) >= minLength) and (len(elem) <= maxlength)):
				seqLen = len(elem)
				if(op == 1):
					if(seqLen not in outDict):
						#A,T,C,G
						tempList = [0,0,0,0]
						outDict[seqLen] = placeNucleotide(tempList, elem)
					elif(seqLen in outDict):
						tempList = outDict[seqLen]
						outDict[seqLen] = placeNucleotide(tempList, elem)
				elif(op == 2):
					placeNucleotide(outList, elem)

	if(op == 1):
		return outDict
	elif(op == 2):
		return outList

#getNucleotideDistro(inFile, 2, 18, 28)
if(op == 1):
	oDict = getNucleotidePosition(inFile, op, minLen, maxLen)

	nucleotideDict = {0:"A", 1:"U", 2:"C", 3:"G"}

	outHandle.write("Length\tA\tU\tC\tG\n")
	outHandleP.write("Length\tA\tU\tC\tG\n")
	outHandlePlot.write("length\tnucleotide\tabundance\n")
	for k,v in sorted(oDict.iteritems()):
               	outHandleP.write(str(k) + "\t")
		outHandle.write(str(k) + "\t")
               	for j,item in enumerate(v):
			percentValue = 0.00
			try:
				percentValue = float(item)/(sum(v)) * 100
			except:
				percentValue = 0.00
                	outHandleP.write(str(percentValue) + "\t")
			outHandle.write(str(item) + "\t")
			nucleotide = nucleotideDict[j]
			outHandlePlot.write(str(k) + "\t" + str(nucleotide) + "\t" + str(percentValue) + "\n")
		outHandleP.write("\n")
		outHandle.write("\n")

elif(op == 2):
	outList = getNucleotidePosition(inFile, op, minLen, maxLen)
	outHandle.write("A\tU\tC\tG\n")
	for item in outList:
		outHandle.write(str(item) + "\t")
