#!/bin/bash
#import the config file here too
dir=`pwd`
realpath=`readlink -f "$0"`
base=${realpath%/*}
dsSeqBase=${base%/*}
root=$dsSeqBase"/pipelineScripts/"

source $dsSeqBase"/pipeline/configFile.txt"
dataPath=$dsSeqBase"/pipeline/data"
function createLogFile()
{
	#creates a log file using the suplied name
	touch $1
}

function logAction()
{
	#Write to a log file different steps of the analysis process: inputs are $1 logfilename and $2 text to write

	echo "$2" >> $1
}

function tick()
{
	#create checkpoint file
	touch $1
}

function checkCommandSuccess()
{
	#Not yet fiunctional...keep in the cooler
	if [[ $? -eq 1 ]]; then
    		echo "Command: ".$1." failed"
	else
		echo "Command: ".$1." succeeded"
	fi
}

function checkFileTypeSpecific()
{
        ftype=`file --mime-type $1 | grep -w "$2" | wc -l`

        if [[ $ftype -gt 0 ]];
        then
                #ensure the file is not empty
                numline=`zcat $1 | head -n 4 | wc -l`
                if [[ $numline -lt 4 ]];
                then
                        #echo file format not correct
                        echo "incorrect fastq format"
                else
                        #check that it is a fastq file
                        fline=`zcat $1 | head -n 1`
                        if [[ ${fline:0:1} == "$3" ]];
                        then
                                #echo "Yes: Gzip & fastq";
                                echo "11"
                        else
                                echo "10"
                                #echo "No: Gzip & Not a fastq file";
                        fi
                fi
        else
                #echo "no"
                fline=`head -n 1 $1`
                if [[ ${fline:0:1} == "$3" ]];
                then
                        echo "01"
                        #echo "yes: Not gzip but fastq";
                else
                        echo "00"
                        #echo "No: Not gzip Not a fastq file";
                fi
        fi
}

function prepRefSeq()
{
	#This fucntion is to prepare the reference sequence usiing Bowtie
	#takes in the fasta sequence of the reference and the name of the index
	#usage will be prepRefSeq "fastaSeq" "IndexName"

	if [[ $mappingSoftware == "Bowtie" ]];
	then
		bowtie-build -f $1 $2
	fi
}

function BowtieMap()
{
	#Mapping a sequence to the reference using Bowtie. It does not return a value as new files will be generated
	#The needed parameters include the number of mismatches $1, the reference index $2, the query sequence file $3 and the output fileName with location $4
	#BowtieMap "mismatches" "refIndex" "InputSeqFile" "OutputFileName.sam"
	#bowtie -a --best -v $1 $2 -f $3 -S -p 12 > $4 2>> $5 &
	bowtie -a --best -v $1 $2 -f $3 -S -p $numCores > $4 2>> $5
	#`disown`
	#bowtie -a --best -v 2 $dir/$4/ref -q $dir/$4/$4_18-28.fastq -S -p 6 > $dir/$4/mapping_$4.sam
}

function mirDeepP()
{
	#$2 is ref genome index | $f is query file | $fname is the name to label the output with (based on the sample ID) | $4 is annotation file gff | $5 is the chromosome length file
        #$1 is the query file | $2 is the reference index | $3 is the reference genome fasta sequence | $4 is annotation file | $5 chromosome length file | $6 is the output name | $7 family size | $8 screen output
	#$9 is an output folder to copy predicted sequences to | $10 is the number of flanking nucleotides to excise
	mirDeepPDir=$root"/miRDP1.3"

	bowtie -a -v 0 "${2}" -f "${1}" > "${6}".aln 2>> "${8}"
	perl $mirDeepPDir/convert_bowtie_to_blast.pl "${6}".aln "${1}" "${3}" > "${6}".bst 2>> "${8}"
	perl $mirDeepPDir/filter_alignments.pl "${6}".bst -c "${7}" > "${6}".filter35.bst 2>> "${8}"
	perl $mirDeepPDir/overlap.pl "${6}".filter35.bst "${4}" -b > "${6}".filter35_overlap 2>> "${8}"
	perl $mirDeepPDir/alignedselected.pl "${6}".filter35.bst -g "${6}".filter35_overlap > "${6}".filter35_overlap.bst 2>> "${8}"
	perl $mirDeepPDir/filter_alignments.pl "${6}".filter35_overlap.bst -b "${1}" > "${6}".filtered35.fa 2>> "${8}"
	perl $mirDeepPDir/excise_candidate.pl "${3}" "${6}".filter35_overlap.bst "${10}" > "${6}".precursors.fa 2>> "${8}"
	cat "${6}".precursors.fa | RNAfold --noPS > "${6}".structures 2>> "${8}"
	bowtie-build -f "${6}".precursors.fa "${6}".precursors 2>> "${8}"
	bowtie -a -v 0 "${6}".precursors -f "${6}".filtered35.fa > "${6}".precursors.aln 2>> "${8}"
	perl $mirDeepPDir/convert_bowtie_to_blast.pl "${6}".precursors.aln "${6}".filtered35.fa "${6}".precursors.fa > "${6}".precursors.bst 2>> "${8}"
	sort +3 -25 "${6}".precursors.bst > "${6}".signatures 2>> "${8}"
	perl $mirDeepPDir/miRDP.pl "${6}".signatures "${6}".structures > "${6}".predictions 2>> "${8}"
	perl $mirDeepPDir/rm_redundant_meet_plant.pl "${5}" "${6}".precursors.fa "${6}".predictions "${6}".nr_predictions "${6}".filter_P_prediction 2>> "${8}"
	cp "${6}".filter_P_prediction "${9}"
}


function checkIndex()
{
        num=`ls $1 | grep $2 | wc -l`

        if [[ $num -gt 0 ]]
        then
                echo "yes"
        else
                echo "no"
        fi
}

function checkCheckPoint()
{
        num=`ls $1 | grep $2 | wc -l`

        if [[ $num -gt 0 ]]
        then
                echo "yes"
        else
                echo "no"
        fi
}

function validateFileReplicates()
{
	#to check that the number of files in here are same as the number of replicates specified in the config file
	sumRep=``
}

function getArrayCombination()
{
	#get the combination of samples to use for the matrix calculations
	#It takes in an array containing the sample label, makes a unique array and then does all the possible combinations
	#Calling the function use the following | testAr=`getArrayCombination "${testArray[@]}"`
	declare -a outArray
        inArray=("$@")
        uniq=($(printf '%s\n' "${inArray[@]}" | sort -u))

        for i in "${!uniq[@]}";
        do
                for j in "${!uniq[@]}";
                do
                        if [ $j -gt $i ];
                        then
                                outArray+=(${uniq[i]}_${uniq[j]})
                        fi
                done
        done
        echo ${outArray[@]}

}

#Generate the matrix used as input into edgeR
function generateMatrix()
{
	#There may need to place another validation of the existence of the file here before execution. What if the file has been moved since the start of the container before getting here???
	#Inputs will be the conditions to be examined - $1 condition 1 and $2 condition 2 | numRep is the number of replicates for each sample
	#check if the number of replicates is 2 or 3

	if [ $3 == 3 ];
	then
		#$dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq"
		#create the different expression headers to be needed here matrixFile=Pipeline/step3_workflow_3replicates_seq.pl
		perl $root/"step3.pl" $1"_1_clip_all"$min_length"-"$max_length".txt" $1"_2_clip_all"$min_length"-"$max_length".txt" $1"_3_clip_all"$min_length"-"$max_length".txt" $2"_1_clip_all"$min_length"-"$max_length".txt" $2"_2_clip_all"$min_length"-"$max_length".txt" $2"_3_clip_all"$min_length"-"$max_length".txt" $1"_"$2"_expMatrix.txt"
		#perl $matrixFile $1"_1_clip_all18-28.txt" $1"_2_clip_all18-28.txt" $1"_3_clip_all18-28.txt" $2"_1_clip_all18-28.txt" $2"_3_clip_all18-28.txt" $2"_5_clip_all18-28.txt" $1"_"$2".expMatrix.txt"
	elif [ $3 == 2 ];
	then
		#create the different expression headers to be needed here
		#step3_2_rep
		perl $root/"step3_2_rep.pl" $1"_1_clip_all"$min_length"-"$max_length".txt" $1"_2_clip_all"$min_length"-"$max_length".txt" $2"_1_clip_all"$min_length"-"$max_length".txt" $2"_2_clip_all"$min_length"-"$max_length".txt" $1"_"$2"_expMatrix.txt"
	elif [ $3 == 1 ];
	then
		#generate keyFile by combining all the sequences from all samples and supply as parameter 4
		python $root/"matrixNoRep.py" 2 $4 $1"_"$2"_expMatrix.txt" $1"_1_clip_all"$min_length"-"$max_length".txt" $2"_1_clip_all"$min_length"-"$max_length".txt"
	fi



	#generate the file without the header
	#`tail -n +2 $1"_"$2".expMatrix.txt" > $1"_"$2".expMatrix_nH.txt"`
	#remove the quotes around the sequence
	#`sed 's/\"//g' $1"_"$2".expMatrix_nH.txt" > $1"_"$2".expMatrix_nHQ.txt"`
	#On completion, generate the checkpoint file
}

function rExpression()
{
	#Run the expression analysis here using the matrix file generated from above
	#$1 the matrix file | $2 sample1 name | $3 sample 2 name | $4 output file | $5 num replicate sample1 | $6 num replicate sample2 | $7 p-value cut-off
	#Rscript $edgeRFile $matrixFile "samp1" "samp2" $outfile $numRep1 $numRep2 $pValEdgeR
	`Rscript $root/"edgeR.r" $1 $2 $3 $4 $5 $6 $7`
}


#########################	FUNCTIONS FOR VALIDATION	#####################################
function fileExist()
{
	#simply check the existence of a file passed as a parameter
	if [ -f "$1" ];
	then
		echo "yes"
	else
		echo "no"
	fi
}


function fileReadable()
{
	#simply check the readability of a file passed as a parameter
        if [ -r "$1" ];
        then
                echo "yes"
        else
                echo "no"
        fi

}

function fileNonEmpty()
{
        #simply check that the file has non-zero size
        if [ -s "$1" ];
        then
                echo "yes"
        else
                echo "no"
        fi
}

function isDirectory()
{
	#simply check that the directory exist
        if [ -d "$1" ];
        then
                echo "yes"
        else
                echo "no"
        fi

}

function checkFileExtension()
{
	#check the extension $2 of a file $1
	#IFS=. read fileName extension <<< $1  #this is the method that does not take care of a file that has names separated by "."

	IFS=. read -r -a fileArray <<< $1
	extension="${fileArray[@]: -1:1}"

	echo $extension;

	if [ $extension = "$2" ];
	then
		echo "yes"
	else
		echo "no"
	fi
}

function uniValidateFile()
{
	#validate a file by checking all input parameters
	numPar=$#

	if [[ $numPar = 2 ]];
	then
		fe=`fileExist $1`

		if [[ $fe = "no" ]];
		then
			echo 5
		else
			fext=`checkFileExtension $1 $2`
			if [[ $fext = "no" ]];
			then
				echo 6
			else
				fr=`fileReadable $1`
				if [[ $fr = "no" ]];
				then
					echo 7
				else
					fne=`fileNonEmpty $1`
					if [[ $fne = "no" ]];
					then
						echo 8
					else
						echo 10
					fi
				fi
			fi
		fi
	else
		echo "incomplete"
	fi
}

function selectOne()
{
	if [[ $1 == "yes" || $2 == "yes" ]];
	then
		echo "yes"
	else
		echo "no"
	fi
}

function transferFile()
{
	#$1 source $@ destination
	fNonEmpty=`fileNonEmpty $1`

	if [[ $fNonEmpty == "yes" ]];
	then
		mv $1 $2
	else
		rm $1
	fi
}
