#File to extract node pairs [edge] for graph construction
#modification for post-grid analysis

import sys

infile = sys.argv[1]
plant_code = sys.argv[2]
outfile = sys.argv[3]		#seqs mapped to plant in miRBase
outfile2 = sys.argv[4]		#seqs not mapped to plant in miRBase
outfile3 = sys.argv[5]		#statistics on seqs

'''class prof_object:

        def __init__(self,id,start,end,value):
                self.id = id
                self.start = start
                self.end = end
                self.value = value
'''
translation_dict = {}
ID_seq_dict = {}
#read the plant translation file here
with open(plant_code) as pc:
	for lines in pc:
		pc_elems = lines.split("\t")
		code = pc_elems[0]
		sci_name = pc_elems[2]
		translation_dict[code] = sci_name

masterList = {} #create an empty dictionary

outHandle = open(outfile,'w') 
outHandleNP = open(outfile2, 'w')
outHandleStat = open(outfile3, 'w')

with open(infile) as f:
	for eachLine in f:
		elems = eachLine.split("\t")				# split each element 
		seq = elems[0]
		seq_count = elems[1]
		ID = elems[0]
		miRNA_ID = elems[9]
		mirbase_ID = elems[2]
		mirbase_ID_cut = mirbase_ID[:11]
		mirbase_ID_elems = mirbase_ID.split("-")
		org_code = mirbase_ID_elems[0]
		fam_code_raw = mirbase_ID_elems[1]
		if(fam_code_raw == "miR"):
			fam_code_raw = mirbase_ID_elems[1] + mirbase_ID_elems[2]
		fam_code = ''.join(i for i in fam_code_raw if i.isdigit())


		if(miRNA_ID in ID_seq_dict):
                        ID_seq_dict[miRNA_ID].append(seq)
                else:
                        ID_seq_dict[miRNA_ID] = [seq]
	
		'''if(seq in ID_seq_dict):
                        ID_seq_dict[seq].append(miRNA_ID)
                else:
                        ID_seq_dict[seq] = [miRNA_ID]'''

		if(seq not in masterList):
			fam_dict = {}
			org_list = []
			org_list.append(org_code)
			fam_dict[fam_code] = org_list
			masterList[seq] = fam_dict
		elif(seq in masterList):
			temp_fam_dict = masterList[seq]
			if(fam_code in temp_fam_dict):
				temp_org_list = temp_fam_dict[fam_code] 
				temp_org_list.append(org_code)
				fam_dict[fam_code] = org_list
				masterList[seq] = fam_dict
			else:
				org_list = []
                        	org_list.append(org_code)
                        	fam_dict[fam_code] = org_list
                        	masterList[seq] = fam_dict


for key, value in masterList.iteritems():
	outHandle.write(ID_seq_dict[key][0] + "\n")
	outHandleNP.write(ID_seq_dict[key][0] + "\n")
	#outHandle.write(str(key) + "\n")
        #outHandleNP.write(str(key) + "\n")

	#fam_dict = value
	for k,v in value.iteritems():
		outHandle.write("miR" + str(k) + "\t")
		outHandleNP.write("miR" + str(k) + "\t")
		#outHandleStat.write(str(key) + "\t" + )
		v_set = set(v)
		v = list(v_set)
		id = "nil"
		mana = "miR" + str(k)
		if(key in ID_seq_dict):
			id = ID_seq_dict[key][0]
		outHandleStat.write(str(key) + "\t" + mana + "\t" + id + "\t" + str(len(v)) + "\n")

		for elements in v:
			real_Name = "Non-plant"
			if(str(elements) in translation_dict):
				real_Name = translation_dict[str(elements)]
				outHandle.write(real_Name + "\t")
			else:
				outHandleNP.write(real_Name + "\t")

		outHandle.write("\n")
		outHandleNP.write("\n")
	outHandle.write("\n\n")
	outHandleNP.write("\n\n")
