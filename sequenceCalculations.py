#!/bin/python
from __future__ import division
import sys

#Sample  Raw_Reads(Library_size) Total_Cleaned_Reads     Filtered_Redundant_Reads(18-35) Filtered_Non-Redundant_Reads(18-35)     Rfam-aligned    Host-aligned(Redundant) Host-aligned(Non-redundant)     Pathogen-aligned(Redundant)     Pathogen-aligned(Non-redundant) Subject-pathogen-aligned(Redundant)     Subject-pathogen-aligned(Non-redundant) Conserved(All)  Conserved(Host)
#ERM1_1  23176392        22473339        10829160        664338  0       196410  45900   6875765 367726          367726  1       0

inFile = sys.argv[1]
outFile = sys.argv[2]

outHandle = open(outFile, 'w')
outHandle.write("Sample\tPercentageHost-Aligned-Redundant-Lib\tPercentageHost-Aligned-Non-Redundant-Lib\tPercentageHost-Aligned-Redundant\tPercentageHost-Aligned-Non-Redundant\tPercentagePathogen-Aligned-Redundant\tPercentagePathogen-Aligned-Non-Redundant")
outHandle.write("\n")

def calculatePercentage(num1,num2):
	percentRes = 0.0
	#print(str(num1) + " and " + str(num2) + " ...")
	try:
		percentRes = (num1/num2) * 100.00
	except:
		percentRes = 0.0

	return percentRes

with open(inFile) as in_f:
	in_f.next()
	for line in in_f:
		lineElems = line.split("\t")
		sample = lineElems[0]
		librarySize = int(lineElems[1])
		cleanedReads = int(lineElems[2])
		filteredReadsRed = int(lineElems[3])
		filteredReadsNonRed = int(lineElems[4])
		hostAlignedReadsRed = int(lineElems[6])
		hostAlignedReadsNonRed = int(lineElems[7])
		PathogenAlignedReadsRed = int(lineElems[8])
		PathogenAlignedReadsNonRed = int(lineElems[9])


		percentHostAlignedRed = calculatePercentage(hostAlignedReadsRed, filteredReadsRed)
		percentHostAlignedNonRed = calculatePercentage(hostAlignedReadsNonRed, filteredReadsNonRed)

		percentHostAlignedRedLibrary = calculatePercentage(hostAlignedReadsRed, librarySize)
                percentHostAlignedNonRedLibrary = calculatePercentage(hostAlignedReadsNonRed, librarySize)

		percentPathogenAlignedRed = calculatePercentage(PathogenAlignedReadsRed, filteredReadsRed)
		percentPathogenAlignedNonRed = calculatePercentage(PathogenAlignedReadsNonRed, filteredReadsNonRed)

		outHandle.write(sample + "\t" + str(percentHostAlignedRedLibrary) + "\t" + str(percentHostAlignedNonRedLibrary) + "\t" + str(percentHostAlignedRed) + "\t" + str(percentHostAlignedNonRed) + "\t" + str(percentPathogenAlignedRed) + "\t" + str(percentPathogenAlignedNonRed) + "\n")
