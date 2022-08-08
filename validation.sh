#!/bin/bash
#script to validate the files and other input paramters before carrying out the analysis
#modified on May 27 2022 to use base paths
realpath=`readlink -f "$0"`
base=${realpath%/*}
dsSeqBase=${base%/*}
root=$dsSeqBase"/pipelineScripts/"
source $dsSeqBase"/pipeline/configFile.txt"		#for dataPath variable
source $root"/funcs.sh"
source $root"/validateFileProc.sh"
dataPath=$dsSeqBase"/pipeline/data"
#Create the log file first
#create a log file
log=$dsSeqBase"/pipeline/logFile.log"
createLogFile $log

#source configFile.txt		#the config file containing the variables
#source funcs.sh			#the file containin the functions to be used

#validation procedure will entail checking the existence of the stated file, non-emptiness and correct extension - later file will be opened to check first few lines for correct content
#check the specified working directory first

wd=`isDirectory $root`
if [ $wd == "no" ];
then
	logAction "$log" "$sourceDir is not a directory"
	logAction "$log" "Please specify a valid directory in order to proceed..."
else
	#####	SCRIPTS	################################
	##### They must all be 'yes' none of these output variables can be 'nil' nor 'Incorrect number of parameters entered for function! 'uniValidateFileProc''
	##### Cummulate the variable for each section and add all sections together at the end, for better management

	if [ -z "$adapterSequence" ]
        then
                adapterSeq="no"
        else
                adapterSeq="yes"
        fi
        #echo $adapterSeq

        logAction "$log" "Validating presence of adapter sequence ..."
        logAction "$log" "$adapterSeq"
	#referenceGenome | edgeRAnalysis | mirdeepAnalysis | microbeAnalysis | mirBaseAnalysis

		#####   GENOME  ################################
	if [[ $referenceGenome == "yes" ]];
	then
	        ref1=`uniValidateFileProc $dataPath/$referenceGenome "fasta"`
	        ref2=`uniValidateFileProc $dataPath/$referenceGenome "fa"`
	        #echo "$ref1 and $ref2"
	        logAction "$log" "Validating reference genome file: $referenceGenome ..."
	        #logAction "$log" "$ref1 || $ref2"
		ref=`selectOne $ref1 $ref2`
	        #echo $ref
		logAction "$log" "Validating reference genome file: $ref"
	else
		ref="yes"
	fi

		#####	MIRDEEP-P	########################
	if [[ $mirdeepAnalysis == "yes" ]];
	then
		mDpGenRef1=`uniValidateFileProc $dataPath/$mirReferenceSequence "fasta"`
		mDpGenRef2=`uniValidateFileProc $dataPath/$mirReferenceSequence "fa"`
		mDpGenRef=`selectOne $mDpGenRef1 $mDpGenRef2`
		logAction "$log" "Validating miRDeepP reference file: $mirReferenceSequence ..."
	        logAction "$log" "$mDpGenRef1 || $mDpGenRef2"


		#echo $mDpGenRef
		mDpAnotation1=`uniValidateFileProc $dataPath/$mirAnnotationFile "gff"`
		mDpAnotation2=`uniValidateFileProc $dataPath/$mirAnnotationFile "gff3"`
		mDpAnotation=`selectOne $mDpAnotation1 $mDpAnotation2`
		logAction "$log" "Validating miRDeepP annotation file: $mirAnnotationFile ..."
	        logAction "$log" "$mDpAnotation1 || $mDpAnotation2"
		echo $mDpAnotation1
		echo $mDpAnotation2
		echo $mDpAnotation
		mDpchromLen=`uniValidateFileProc $dataPath/$mirChromosomeLength "txt"`
		logAction "$log" "Validating chromosome length file: $mirChromosomeLength ..."
	        logAction "$log" "$mDpchromLen"

		#echo $mDpchromLen
		if [[ $mDpGenRef == "yes" && $mDpAnotation == "yes" && $mDpchromLen == "yes" ]];
	        then
	                mDpfiles="yes"
			logAction "$log" "Successfully validated all mirDeep-P requirements!"
	        else
	                mDpfiles="no"
			logAction "$log" "Failed to successfully validate all mirDeep-P requirements!"
	        fi
		#echo $mDpfiles
	else
		mDpfiles="yes"
	fi

		#####	MIRBASE		########################
	if [[ $ == "yes" ]];
	then
		mirRef1=`uniValidateFileProc $dataPath/$mirBaseReference "fasta"`
	        mirRef2=`uniValidateFileProc $dataPath/$mirBaseReference "fa"`
	        mirRef=`selectOne $mirRef1 $mirRef2`
	        #echo $mirRef
		logAction "$log" "Validating reference sequence file for miRBase: $mirBaseReference ..."
	        logAction "$log" "$mirRef1 || $mirRef2"

		logAction "$log" "Validating number of mismatch for miRBase mapping ..."
	        if [ -z "$mirBaseMismatch" ]
	        then
	                misMir="no"
	        else
	                misMir="yes"
	        fi
		#echo $misMir
		logAction "$log" "$misMir"

		f1NonEmpty=`fileNonEmpty $dataPath/$mirBaseReference`
                if [[ $f1NonEmpty == "yes" ]];
                then
                        `awk '{if(NR%2==0){gsub("U","T",$0); print;}else{print}}' $dataPath/$mirBaseReference > $dataPath/"mirBaseRef.tmp"`
                        f2NonEmpty=`fileNonEmpty $dataPath/"mirBaseRef.tmp"`
                        if [[ $f2NonEmpty == "yes" ]];
                        then
                                #move file for expression analysis
                                mv $dataPath/"mirBaseRef.tmp" $dataPath/$mirBaseReference
                                convertMirRefFile="yes"
                        else
                                convertMirRefFile="no"
                                logAction "$log" "Temporary miRBase reference file empty!"
                        fi
                else
                        logAction "$log" "Original miRBase reference file empty!"
                        convertMirRefFile="no"
                fi


                if [[ $mirRef == "yes" && $misMir == "yes" && $convertMirRefFile == "yes" ]];
                then
                        mirSect="yes"
                        logAction "$log" "Successfully validated all mirBase mapping requirements!"

	        else
	                mirSect="no"
			logAction "$log" "Failure validating all mirBase mapping requirements!"
	        fi

		#echo $mirSect
	else
		mirSect="yes"

	fi

		#####	EDGE-R		########################
		#checking that the p-value is input
	if [[ $edgeRAnalysis == "yes" ]];
	then
		logAction "$log" "Validating pValue setting for edgeR analysis ..."
		if [ -z "$pValEdgeR" ]
	        then
	                edgerPVal="no"
	        else
	                edgerPVal="yes"
	        fi
		#echo $edgerPVal
		logAction "$log" "$edgerPVal"
	else
		edgerPVal="yes"
	fi
		#####	MICROBE		########################
	if [[ $microbeAnalysis == "yes" ]];
	then
		virRef1=`uniValidateFileProc $dataPath/$microbeReference "fasta"`
		virRef2=`uniValidateFileProc $dataPath/$microbeReference "fa"`
		virRef=`selectOne $virRef1 $virRef2`
		#echo $virRef
		logAction "$log" "Validating reference sequence file for microbe repository: $microbeReference ..."
	        logAction "$log" "$virRef1 || $virRef2"

		logAction "$log" "Validating microbe-of-interest setting for microbe analysis ..."
		if [ -z "$subjectMicrobe" ]
	  	then
			subVir="no"
		else
			subVir="yes"
			logAction "$log" "$subjectMicrobe selected as microbe of interest"
		fi
		#echo $subVir

		logAction "$log" "Validating number of mismatch for microbe repository mapping ..."

		if [ -z "$microbeMismatch" ]
	        then
	                misVir="no"
	        else
	                misVir="yes"
	        fi
		#echo $misVir
		logAction "$log" "$misVir"

		if [[ $virRef == "yes" && $misVir == "yes" ]];
	        then
	                virSect="yes"
			logAction "$log" "Successfully validated all microbe repository mapping requirements!"
	        else
	                virSect="no"
			logAction "$log" "Failure validating all microbe repository mapping requirements!"
	        fi
		#echo $virSect
	else
		virSect="yes"
	fi
fi

#echo "adapter seuquence -"$adapterSeq
#echo "host section - "$ref
#echo "microbe section - "$virSect
#echo "edgeR pValue - "$edgerPVal
#echo "mirBase section - "$mirSect
#echo "mirDeep-P section - "$mDpfiles

function checkValidate()
{
	#combine all results
	if [[ $adapterSeq == "yes" && $ref == "yes" && $mDpfiles == "yes" && $mirSect == "yes" && $virSect == "yes" && $edgerPVal == "yes" ]];
	then
		logAction "$log" "All validation successful!"
		echo "yes"
	else
		logAction "$log" "Failed to successfully validate all files!"
		echo "no"
	fi
}
