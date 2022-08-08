#group mirbase sequences into families
#	AUTHOR: OLAGUNJU A. TEMITAYO	OCTOBER 2019
#python A1/format_mirname_api.py E2/E2.mir.txt /mnt/tayo/testSRNAPipeStubDocker/pipelineScripts/viridiplantae_plants.txt E2/E2_mirname_api.txt
import sys

infile = sys.argv[1]
plant_code = sys.argv[2]
outfile = sys.argv[3]		#seqs mapped to plant in miRBase

def hasNumber(s):
    return any(i.isdigit() for i in s)

translation_dict = {}
ID_seq_dict = {}
actual_ID = {}

#read the plant translation file here
with open(plant_code) as pc:
	for lines in pc:
		pc_elems = lines.split("\t")
		code = pc_elems[0]
		sci_name = pc_elems[2]
		translation_dict[code] = sci_name

masterList = {} #create an empty dictionary

outHandle = open(outfile,'w') 

with open(infile) as f:
	for eachLine in f:
		elems = eachLine.split("\t")				# split each element 
		seq = elems[9]
		seq_count = elems[1]
		ID = elems[0]
		miRNA_ID = elems[9]
		mirbase_ID = elems[2]
		
		mirbase_ID_elems = mirbase_ID.split("-")
		org_code = mirbase_ID_elems[0]
		if(len(mirbase_ID_elems) == 2):
			sec = mirbase_ID_elems[1]
			if(sec.isalpha()):
				payload = sec
			elif(hasNumber(sec)):
				if("." in sec):
					sec_elems = sec.split(".")
					payload = ''.join(i for i in sec_elems[0] if i.isdigit())
				elif("." not in sec):
					payload = ''.join(i for i in sec if i.isdigit())
		elif(len(mirbase_ID_elems) == 3):
			sec = mirbase_ID_elems[1]
			if(sec != "miR"):
				if(sec.isalpha() and sec != "bantam"):
					third = mirbase_ID_elems[2]
					num = ''.join(i for i in third if i.isdigit())
					payload = str(sec) + "-" + str(num)
				elif(sec.isalpha() and sec == "bantam"):
					payload = sec
				elif(hasNumber(sec)):
					payload = ''.join(i for i in sec if i.isdigit())
			elif(sec == "miR"):
				third = mirbase_ID_elems[2]
				payload = ''.join(i for i in third if i.isdigit())
		elif(len(mirbase_ID_elems) == 4):
			sec = mirbase_ID_elems[1]
			if(sec != "miR"):
				if(sec.isalpha()):      #(dsi-let-7a) | sme-bantam-c-3p
                                        third = mirbase_ID_elems[2]
                                        num = ''.join(i for i in third if i.isdigit())
					if(len(num) > 0):
                                        	payload = str(sec) + "-" + str(num)
					else:
						payload = str(sec)
				elif(hasNumber(sec)):
					payload = ''.join(i for i in sec if i.isdigit())
			elif(sec == "miR"):
				third = mirbase_ID_elems[2]
				payload = ''.join(i for i in third if i.isdigit())
		org_name = "no record"
		if(org_code in translation_dict):
			org_name = translation_dict[org_code]
				
		#outHandle.write(str(mirbase_ID) + "\t-\t" + str(payload) + "\t" + str(org_name) + "\n")
		
		if(payload not in masterList):
			masterList[payload] = [org_name]
			if(payload not in ID_seq_dict):
				ID_seq_dict[payload] = [elems[9]]
			else:
				ID_seq_dict[payload].append(elems[9])
			#masterList[payload] = [org_name]
		elif(payload in masterList):
			masterList[payload].append(org_name)
			if(payload not in ID_seq_dict):
                                ID_seq_dict[payload] = [elems[9]]
                        else:
                                ID_seq_dict[payload].append(elems[9])


		
for key, value in masterList.iteritems():
	outHandle.write("miR-" + str(key) + "\t" + str(ID_seq_dict[key][0]) + "\t")

	val = set(value)
	value = list(val)
        
	for i,elements in enumerate(value):
		if(i == (len(value) -1)):
			outHandle.write(elements)
		else:
			outHandle.write(elements + " _ ")
	outHandle.write("\n")
