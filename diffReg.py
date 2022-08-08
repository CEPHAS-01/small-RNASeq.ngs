#AUTHOR: OLAGUNJU A. TEMITAYO 10 OCT 2020
import sys

if(len(sys.argv) < 7):
	print("Please use correct number of arguments")
	usageText='''\nUSAGE: python diffReg.py inputFile outputFile Operation logFC pVal logFCColumn pValColumn\n'''
	print(usageText)
	sys.exit()
else:
	in_file = sys.argv[1]
	out_file = sys.argv[2]
	op = sys.argv[3]
	logFC = float(sys.argv[4])
	pVal = float(sys.argv[5])
	logFCCol = int(sys.argv[6]) - 1
	pValCol = int(sys.argv[7]) - 1

	outMapHandle = open(out_file,'w')
		
	with open(in_file) as k:
		for kLines in k:
			keyElems = kLines.split("\t")
			if(len(keyElems) > 7):
				pValue = float(keyElems[pValCol])
				fChange = float(keyElems[logFCCol])
		
				if(op == "1"):
					if(((fChange >= logFC) or (fChange <= -logFC)) and (pValue < pVal)):
						outMapHandle.write(kLines)
				elif(op == "2"):
					if((fChange >= logFC) and (pValue < pVal)):
                                                outMapHandle.write(kLines)
				elif(op == "3"):
					if((fChange <= -logFC) and (pValue < pVal)):
                                                outMapHandle.write(kLines)
	#print("Processing completed!")
	
