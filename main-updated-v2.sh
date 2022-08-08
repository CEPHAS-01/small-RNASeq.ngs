#!/bin/bash
#!/bin/env Rscript
################################################################################################################
#
#	Temitayo A. Olagunju July 2019
#
#	Command line: sh small_RNA_Seq_Pipeline.sh "logFile.log" &
#modified for the new addition to statistics file
# May 27 2022: Modified the host mapping to produce unmapped reads file
#################################################################################################################
start=`date +%H:%M:%S`
START_1=$(date +%s)

#Check and set the source directory
if [ -z "$sourceDir" ]
then
	dir=`pwd`
else
	dir=$sourceDir
fi

realpath=`readlink -f "$0"`
base=${realpath%/*}
dsSeqBase=${base%/*}
root=$dsSeqBase"/pipelineScripts/"

#root="/pipelineScripts"
source $dsSeqBase"/pipeline/configFile.txt"
source $root"/funcs.sh"
source $root"/validation.sh"

dataPath=$dsSeqBase"/pipeline/data"
UserDir=$sourceDir
#dataPath
#create a log file
log=$dir"/logFile.log"
#createLogFile $log	#this has been created in validation.sh file, the ID only being used for reference here

#validate here and decide whether to continue the analysis or not
chkValidate=`checkValidate`

if [[ $chkValidate == "yes" ]];
then
	echo "Script execution started at time $start... == $START_1" >> $log

	#create file for mapping stderr output to screen
	mappinglog=$dir"/mappingLogFile.log"
	createLogFile $mappinglog

	#validate here and decide whether to continue the analysis or not
	chkValidate=`checkValidate`

	#create a directory to hold the sequence segregation files
	seqSegregate=$dir"/seqSeg"
	mkdir -p $seqSegregate

	#create a file to hold the header of the edgeR exression output
	expHeader=$seqSegregate"/expHeader.txt"
	createLogFile $expHeader

	#Create a directory to hold the checkpoint files
	checkPoint=$dir"/checkpoint"
	mkdir -p $checkPoint

	#Create a directory to keep mapping indices
	indexDir=$dir"/Index"
	mkdir -p $indexDir

	#create a permanent directory for each data to be retained in the case of a cleanup of intermediate files when the analysis completes
	permDir=$dir"/sRNAOutput"
	mkdir -p $permDir

	#create directory for host-mapped sequences
	hostMapDir=$dir"/hostMap"
	mkdir -p $hostMapDir

	#create a directory for the plots
	plotsDir=$permDir"/Plots"
	mkdir -p $plotsDir

	#create a directory for the reports
	reportsDir=$permDir"/Reports"
	mkdir -p $reportsDir

	#create directory for the mapping files
	mappingFilesDir=$permDir"/AlignmentMaps"
	mkdir -p $mappingFilesDir

	refBuilt="no"

	#ceate a directory for mirdeep-P analysis
	if [[ $mirdeepAnalysis == "yes" ]];
	then
		mirDeepP=$dir"/mirDeepP"
		mkdir -p $mirDeepP

		mirDeepP_2=$dir"/mirDeepP/mDP2"
		mkdir -p $mirDeepP_2

		mirDeepPOutput=$permDir"/NovelPrediction"
		mkdir -p $mirDeepPOutput
	fi

	#create a directory for the mirbase family grouping
	if [[ $mirBaseAnalysis == "yes" ]];
        then
                mirBaseFamily=$dir"/mirFamily"
                mkdir -p $mirBaseFamily

		mirBaseOutput=$permDir"/miRBase"
		mkdir -p $mirBaseOutput
        fi


	#Determine whether to carry out expression analysis so as to create a folder to put the files to be used
	if [[ $edgeRAnalysis == "yes" ]];
	then
		#create the directory here and call it expression
		expDirName=$dir"/expression"
		mkdir -p $expDirName
		#get the list of all the samples being used...needed for now in the edgeR matrix expression calculations. This may need to be moved into the conditional clause below to be used only if expression analysis is needed
        	declare -a sampleArray

		expOutput=$permDir"/ExpressionProfile"
		mkdir -p $expOutput

		expReport=$dir/"expressionStatistics.txt"
		rm -f $expReport
	        createLogFile $expReport
        	echo -e "DIFFERENTIAL REGULATION PROFILE\n" >> $expReport
        	echo -e "Samples\tTotal_reported\tDiff_Regulated\tUp-regulated\tDown-regulated" >> $expReport
	fi

	#create a directory for the virus output in the permanent directory if virus analysis is to be done
        if [[ $virusAnalysis == "yes" ]];
        then
                virusOutput=$permDir"/Virus"
                mkdir -p $virusOutput

		virReport=$dir/"virusStats.txt"
		rm -f $virReport
        	createLogFile $virReport
        	echo -e "VIRUS ANALYSIS REPORT\n" >> $virReport
        	echo -e "Sample\tRedundant_Sequence_Count\tNon-redundant_Sequence_Count\tNum_DB_Accessions\tNon-redundant_${subjectVirus}_Count\tRedundant_${subjectVirus}_Count" >> $virReport

        fi

	#CREATE statistics files here
	smallRNADistro=$dir"/smallRNALengthDistribution.txt"
	rm -f $smallRNADistro
	createLogFile $smallRNADistro

	#rawSmallRNADistro=$dir"/rawSmallRNALengthDistribution.txt"
	#rm -f $rawSmallRNADistro
	#createLogFile $rawSmallRNADistro

	echo -e "Length\tCount\tSample" >> $smallRNADistro
	#echo -e "Length\tCount\tSample" >> $rawSmallRNADistro

	readsStat=$dir"/readsStatistics.txt"
	rm -f $readsStat
	createLogFile $readsStat
	statHeader="Sample\tRaw_Reads(Library_size)\tTotal_Cleaned_Reads\tFiltered_Redundant_Reads("$min_length"-"$max_length")\tFiltered_Non-Redundant_Reads("$min_length"-"$max_length")\tRfam-aligned\tHost-aligned(Redundant)\tHost-aligned(Non-redundant)\tPathogen-aligned(Redundant)\tPathogen-aligned(Non-redundant)\tSubject-pathogen-aligned(Redundant)\tSubject-pathogen-aligned(Non-redundant)\tConserved(All)\tConserved(Host)" 
	echo -e $statHeader >> $readsStat

	genStat=$dir"/analysisStatistics.txt"
	rm -f $genStat
	createLogFile $genStat
	echo -e $statHeader >> $genStat

	fileMap=$dir"/sampleFilesMap.txt"
        rm -f $fileMap
        createLogFile $fileMap
	echo -e "ANALYSIS FILE MAP\nOriginal_File_Name\tAnalysis_File_Name\n"

	#The loop for each file in the data folder begins here
	#numFiles=`ls -l $dir/*.fastq | wc -l`
	logAction "$log" "GENERAL:: The number of samples seen in the directory is : $numFiles"

	#move into data directory to create the filemap associative array first
	cd $dataPath
	dataDir=`pwd`

	declare -A FILEMAP

	dirs=`ls -F $dataDir | grep \/$`
	for folder in $dirs
        do
                #$folder # the folder with trailing slash
                fname=${folder%/}       #folder name alone without any slash
                cd $fname
                numFiles=`ls * | grep '.fastq\|.gz'`                     #get a  list of the files with the specified extension

                fileArr=($numFiles)                             #convert the list into an array
                for i in "${!fileArr[@]}";
                do
                        sampNum=`expr $i + 1`
                        fullFileName=${fileArr[i]}                      #filename with extesion
                        fileName=${fullFileName%%.*}                    #filename without extension

                        tmp=$fname"_"$sampNum
                        #echo -e "$file\t$fileName\t$analysisSampleName"
			echo -e "$fullFileName\t$tmp" >> $fileMap
                        FILEMAP[$tmp]=$fileName

			if [[ $edgeRAnalysis == "yes" ]];
                        then
                                sampleArray+=($fname)
                        fi

			logAction "$log" "Ds-Seq: Starting analysis of $tmp..."

                        #check if the file is a direct fastq file #"01"
                        #file | extension pattern (gzip) | file signature (>, @)
                        chkResFQ=`checkFileTypeSpecific $dataDir/$fname/$fullFileName "gzip" "@"`
                        echo $chkResFQ
                        if [[ $chkResFQ == 01 ]];
                        then
                                #use the file directly
                                echo "$file direct fastq file"
				`cp $dataDir/$fname/$fullFileName $dataDir/"tmpInput.fastq"`
                        else
                                #check if the file is a gzip fastq file #"11"
                                chkResGZ=`checkFileTypeSpecific $dataDir/$fname/$fullFileName "gzip" "@"`
                                echo $chkResGZ
                                if [[ $chkResGZ == 11 ]];
                                then
                                        #zcat the file
                                        echo "$file not a direct fastq but a zipped fastq"
					`zcat $dataDir/$fname/$fullFileName > $dataDir/"tmpInput.fastq"`
                                else
                                        #log that not sure what type of file it is and cannot be used...exit the script or ignore this stage
                                        #echo "$file not a fastq and not a zipped fastq"
					echo "File type undefined"
                                fi
                        fi

			#create and use the file here
			file=$dataDir/"tmpInput.fastq"

			#make the directory for the current file if it does not already exist else leave it
			mkdir -p $dir/$tmp

			realSampleName=${FILEMAP[$tmp]}

			#make a length distribution of the raw fastq file prior to commencing analysis
			#awk 'NR%4==2{lengths[length($0)]++} END {for(l in lengths) {print l, lengths[l]}}' $file | sort -k 1 -n | awk '{ print $0 "\t" "'"$realSampleName"'"}' >> $rawSmallRNADistro

			logAction "$log" "READS FILTERING:: step1..."
			perl $root/"step1.pl" $file $adapterSequence $dir/$tmp $tmp 2>> $mappinglog

			#check if the HD protocol is to be used here, if yes, then use it if not, go on to the next step
			if [[ $HDprotocol == "no" ]];
			then
				logAction "$log" "READS FILTERING:: Not using the HD protocol, skipping this step..."
				cp $dir/$tmp/$tmp.fastq $dir/$tmp/$tmp"_clip.fastq"
			else
				logAction "$log" "READS FILTERING:: HD protocol clipping..."
				perl $root/"HD.pl" $dir/$tmp/$tmp.fastq 4 $dir/$tmp/$tmp"_clip.fastq" 2>> $mappinglog
			fi

			#obtain length distribution after adaptor removal
			logAction "$log" "READS FILTERING:: Creating length distribution file after adaptor removal..."
			awk 'NR%4==2{lengths[length($0)]++} END {for(l in lengths) {print l, lengths[l]}}' $dir/$tmp/$tmp"_clip.fastq" | sort -k 1 -n > $dir/$tmp/$tmp"_clip_no-adaptor_length_plot.txt"
			#PENDING::::Make the plot of the length distribution - perhaps to transfer the files first



			logAction "$log" "READS FILTERING:: Step2..."
			perl $root/"step2.pl" $dir/$tmp/$tmp"_clip.fastq" $dir/$tmp $min_length $max_length $freq_cutoff 2>> $mappinglog

			#Reads statistics will begin  here
			smplID=$tmp
			`awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l "\t" lengths[l]}}' $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq" > $dir/$tmp/$tmp"smallRNADistro_raw.txt"`
			#realSampleName=${FILEMAP[$tmp]}
			sort -k1,1 $dir/$tmp/$tmp"smallRNADistro_raw.txt" | awk '{ print $0 "\t" "'"$realSampleName"'"}' >> $smallRNADistro

			#total raw reads count
			#check file exist and non-empty
			grepKey1="total sequences:"
			grepKey2="exp:"

			f0NonEmpty=`fileNonEmpty $dir/$tmp/"log-seq_clean.txt"`
                	if [[ $f0NonEmpty == "yes" ]];
                	then
				#rawReadsLib=`head -n 2 $dir/$tmp/"log-seq_clean.txt" | tail -n 1 | awk -F'[^0-9]*' '$0=$2'`
				rawReadsLib=`grep "${grepKey1}" $dir/$tmp/"log-seq_clean.txt" | awk -F'[^0-9]*' '$0=$2'`
				cleanedReads=`grep "${grepKey2}" $dir/$tmp/"log-seq_clean.txt" | awk -F'[^0-9]*' '$0=$2'`
			else
				cleanedReads="nil"
				rawReadsLib="nil"
			fi

			#redundant reads count
			f1NonEmpty=`fileNonEmpty $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq"`
                	if [[ $f1NonEmpty == "yes" ]];
                	then
				redRawReadsNum=`awk '{print $0}' $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq" | wc -l`
				redReads=`expr $redRawReadsNum / 4`

				awk '(NR%4==2)' $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq" > $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq"

				#create file for matrix generation for no replicates samples
				if [[ $samplesReplicated == "no" ]];
                        	then
					`awk '(NR%4==2)' $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq" > $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq.seq"`
					matNoRepNonEmpty=`fileNonEmpty $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq.seq"`
	                        	if [[ $matNoRepNonEmpty == "yes" ]];
	                        	then
						#move file for expression analysis
						cp $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq.seq" $expDirName
					fi
				fi

			else
				redReads="nil"
			fi

			#non-redundant reads count
			#check file exists and non-empty
			f2NonEmpty=`fileNonEmpty $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".fasta"`
                	if [[ $f2NonEmpty == "yes" ]];
                	then
				nonRedRawReadsNum=`awk '{print $0}' $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".fasta" | wc -l`
				nonRedReads=`expr $nonRedRawReadsNum / 2`
			else
				nonRedReads="nil"
			fi

			RfamMap="nil"
			if [[ $mapToRfam == "yes" ]];
                        then
                                #get the file
                                #prepare Rfam reference
                                referenceRfam=$dataPath"/"$referenceRfam
                                refRfamName=${referenceRfam##*/}
                                refRfamIndex=${refRfamName%%.*}

                                checkRfamIndexVar=`checkIndex $indexDir $refRfamIndex`

                                if [[ $checkRfamIndexVar == "no" ]];
                                then
                                        logAction "$log" "RFAM MAPPING:: Rfam reference Index does not exist...preparing the reference index..."
                                        #mod: instead of moving into the index directory, why not use the absolute path?
                                        cd $indexDir
                                        prepRefSeq $referenceRfam $refRfamIndex
                                        cd $dataDir/$fname
				fi


				#map with bowtie
				logAction "$log" "RFAM MAPPING:: Mapping $tmp to the reference genome..."

        			bowtie -a --best -v $RfamMismatch $indexDir/$refRfamIndex -f $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" --un $dir/$tmp/$tmp"_Rfam_unaligned.fasta" -S -p $numCores > $dir/$tmp/$tmp"_Rfam.sam" 2>> $mappinglog

				mapFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp"_Rfam.sam"`
                                if [[ $mapFileNonEmpty == "yes" ]];
                                then
					#get only the mapped reads
					grep -v @ $dir/$tmp/$tmp"_Rfam.sam" | awk '{if($2 == 0 || $2 == 16) {print $0}}' > $dir/$tmp/$tmp"_Rfam_mapped.txt"

					#get the number of reads mapped to Rfam
					RfamMap=`cut -f 10 $dir/$tmp/$tmp"_Rfam_mapped.txt" | sort | uniq -c | wc -l`
				else
					logAction "$log" "RFAM MAPPING:: Rfam alignment file empty..."
					RfamMap="Empty"
				fi

				#straighten the fasta file to have access to the sequence for the unaligned sequences fastq files
				awk '{if(NR%2 == 0) {print $0}}' $dir/$tmp/$tmp"_Rfam_unaligned.fasta" > $dir/$tmp/$tmp"_Rfam_unaligned_seq.fasta"

				#do this step only if the analysis requires removal of Rfam-mapped reads.
				#change the name of the unaligned fasta file to the one that is to be used for all alignments
				#change the original file so it is not overwritten
				if [[ $excludeRfam == "yes" ]];
                        	then
					logAction "$log" "RFAM MAPPING:: Pipeline excluding Rfam sequences..."
					mv $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length"_old.fasta"

					mv $dir/$tmp/$tmp"_Rfam_unaligned.fasta" $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta"
				else
					logAction "$log" "RFAM MAPPING:: Not excluding Rfam sequences..."
				fi
				#delete the Rfam mapping SAM file
				logAction "$log" "RFAM MAPPING:: Cleaning up Rfam SAM and mapping files for $tmp..."

				#if [[ $keepMapFiles == "yes" ]];
				#then
					#mv $dir/$tmp/$tmp"_Rfam.sam" $mappingFilesDir
				#else
					#rm $dir/$tmp/$tmp"_Rfam.sam"
				#fi
                        fi

			#CREATE FILE FOR NUCLEOTIDE DISTRIBUTION AND PLOT
			awk '{if(NR%2 == 0) {print $0}}' $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" > $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".Plotseq"
			#python nucleotideDistribution.py infile outfile operation(1-dict, 2-list) minLen maxLen
			plotOutFile=$tmp"-nucleotideDistro.txt"
			python $root/"nucleotideDistribution.py" $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".Plotseq" $plotOutFile 1 $min_length $max_length

			#MAKE THE PLOT OF THE OUTPUT
			dqNonEmpty=`fileNonEmpty "percentage-Plot-"$plotOutFile`

        		if [[ $dqNonEmpty == "yes" ]];
        		then
				plotTitle=$tmp" Nucleotides Distribution"
				if [[ $plotFormat == "png" ]];
                                then
					outNucleotidePlot=$dir/$tmp/$tmp"_Nucleotide_Distro_plot.png"
				elif [[ $plotFormat == "pdf" ]];
                                then
					outNucleotidePlot=$dir/$tmp/$tmp"_Nucleotide_Distro_plot.pdf"
				fi
				#outNucleotidePlot=$dir/$tmp/$tmp"_Nucleotide_Distro_plot"
				Rscript $root/"plotNucleotidesDistro.R" "percentage-Plot-"$plotOutFile $outNucleotidePlot $plotFormat $plotTitle

				drNonEmpty=`fileNonEmpty $outNucleotidePlot`
				if [[ $drNonEmpty == "yes" ]];
                        	then
					logAction "$log" "PLOT - Nucleotide distro:: Plot created for $tmp"
					#move the plot to the output folder
					mv $outNucleotidePlot $plotsDir
				fi
				mv "percentage-Plot-"$plotOutFile "percentage-"$plotOutFile $plotOutFile $reportsDir
			else
				logAction "$log" "PLOT - Nucleotide distro:: File for nucleotide distro empty"
			fi

			#move file for expression analysis
                        if [[ $edgeRAnalysis == "yes" ]];
                        then
				#convert the clip_all.txt file from underscore- to tab-delimited
				#exclude the Rfam-mapped reads from the file before moving for expression analysis
				#if in the index exclude, otherwise print it
				if [[ $excludeRfam == "yes" ]];
                                then
					RfamMappedNotEmpty=`fileNonEmpty $dir/$tmp/$tmp"_Rfam_mapped.txt"`

					if [[ $RfamMappedNotEmpty == "yes" ]];
                        		then
						logAction "$log" "RFAM MAPPING:: Processing file for expression analysis export..."
						awk '{gsub("_","\t",$0); print;}' $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".txt" > $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_tabbed.txt"
						awk 'NR==FNR{a[$10]=$0;next} {if($1 in a){next} else {print $0}}' $dir/$tmp/$tmp"_Rfam_mapped.txt" $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_tabbed.txt" > $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_tabbed_less_Rfam.txt"
						#convert back to underscore-delimited
						awk '{gsub("\t","_",$0); print;}' $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_tabbed_less_Rfam.txt" > $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_underscore_less_Rfam.txt"

						mv $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".txt" $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_old.txt"
                                		#move the file for the matrix creation here | using copy for the while
						mv $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length"_underscore_less_Rfam.txt" $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".txt"
					fi
				fi

				logAction "$log" "POST-RFAM MAPPING:: Moving files into expression analysis folder..."
                                cp $dir/$tmp/$tmp"_clip_all"$min_length"-"$max_length".txt" $expDirName

				#cleanup
				#rm $dir/$tmp/$tmp"_Rfam_mapped.txt"
                        fi


			if [[ $mapToReference == "yes" ]];
			then
				genMap="nil"
				genMapRed="nil"
				#Check if the index exists already if it does not, then prepare it, else use it
				referenceGenome=$dataPath"/"$referenceGenome
				refName=${referenceGenome##*/}
        			refIndex=${refName%%.*}

				checkIndexVar=`checkIndex $indexDir $refIndex`

				if [[ $checkIndexVar == "no" ]];
				then
					logAction "$log" "HOST MAPPING:: Genome reference Index does not exist...preparing the reference index..."
					#mod: instead of moving into the index directory, why not use the absolute path?
					cd $indexDir
					prepRefSeq $referenceGenome $refIndex
					refBuilt="yes"
					cd $dataDir/$fname
					#cd ../
				else
					refBuilt="yes"
				fi

				if [[ $refBuilt == "yes" ]];
                                then
					logAction "$log" "HOST MAPPING:: Mapping $tmp to the reference genome..."
					#--un $dir/$tmp/$tmp"_Rfam_unaligned.fasta"
					bowtie -a --best -v $HostMismatch $indexDir/$refIndex -f $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" --un $dir/$tmp/$tmp"_unaligned.fasta" -S -p $numCores > $dir/$tmp/$tmp.sam 2>> $mappinglog
					#copy the unaligned reads file
					cp $dir/$tmp/$tmp"_unaligned.fasta" $permDir
					#BowtieMap $HostMismatch $indexDir/$refIndex $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" $dir/$tmp/$tmp.sam $mappinglog
				else
					logAction "$log" "HOST MAPPING:: Genome reference Index still not found..."
				fi

				#keepMapFiles
				if [[ $keepHostMapFiles == "yes" ]];
	                        then
        	                        #move the file to where it will be processed into BAM files
					#create the directory for host map in the maping folder
					hostMapFolder=$mappingFilesDir"/host"
					mkdir -p $hostMapFolder
					#$fileName			| the actual file name, may not be necessary to use it yet if the files will still be zipped later after all analysis
					#mod: copy only if it is not empty
					mapFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp.sam`
                        		if [[ $mapFileNonEmpty == "yes" ]];
                        		then
                	                	cp $dir/$tmp/$tmp.sam $hostMapFolder
					fi
                        	fi

				#post-mapping processing of files will begin here
				#get only mapped output lines
				#mod: this block should be moved into the part above to be executed only if the file is non-empty
				mapFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp.sam`
                                if [[ $mapFileNonEmpty == "yes" ]];
                                then
					#move the SAM mapping file
					#cp $dir/$tmp/$tmp.sam $hostMapFolder

					grep -v @ $dir/$tmp/$tmp.sam > $dir/$tmp/$tmp.txt
					awk '{if($2 == 0 || $2 == 16) {print $0}}' $dir/$tmp/$tmp.txt > $dir/$tmp/$tmp"_mapped.txt"

					#statistics file recording
                                        genMap=`cut -f 1 $dir/$tmp/$tmp"_mapped.txt" | sort | uniq -c | wc -l`
					#genMapRed=`wc -l $dir/$tmp/$tmp"_mapped.txt" | awk -F'[^0-9]*' '$0=$2'`

					#mod: copy only if non-empty
					#cp $dir/$tmp/$tmp"_mapped.txt" $hostMapDir
					cp $dir/$tmp/$tmp"_mapped.txt" $hostMapDir

					#cut only the sequence and then copy the file to the sequence segregation folder
					`cut -f 10 $dir/$tmp/$tmp"_mapped.txt" | awk '!seen[$1]++' > $dir/$tmp/$tmp"_map_cut.txt"`
					#mod: ensure that all files to be copied or moved are not empty
					#straighten the fastq file for the purpose of counting the abundance of mapped reads
					awk '{if(NR%2 == 0) {print $0}}' $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length".fastq" | awk '{if(NR%2 != 0) {print $0}}' > $dir/$tmp/$tmp"_fastq_"$min_length"-"$max_length".seq"
					#$dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq"
					genMapRed=`awk 'NR==FNR{a[$1]=$0;next} ($1 in a) {print $0}' $dir/$tmp/$tmp"_map_cut.txt" $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq" | wc -l`
					#genMapRed=`awk 'NR==FNR{a[$1]=$0;next} ($1 in a) {print $0}' $dir/$tmp/$tmp"_map_cut.txt" $dir/$tmp/$tmp"_fastq_"$min_length"-"$max_length".seq" | wc -l`
					cp $dir/$tmp/$tmp"_map_cut.txt" $seqSegregate

					#length distribution of reads mapped to the reference host genome with reference to the cleaned reads
					awk 'NR==FNR{a[$0]=$0;next} $0 in a {print $0}' $dir/$tmp/$tmp"_map_cut.txt" $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq" | awk '{lengths[length($0)]++} END {for(l in lengths) {print l, lengths[l]}}' | sort -k 1 -n > $dir/$tmp/$tmp"_post-host-mapping-length-distro-red.txt"
					#PENDING:::: PLOT the distribution
				else
					logAction "$log" "HOST MAPPING:: Empty host mapping SAM file for $tmp ..."
					genMap="Empty"
					genMapRed="Empty"
				fi
				#tick that the hostmap stage for this has been completed
				tick $checkPoint/$tmp"_hostMap.chkpt"
			else
				genMap="No_Map_Option"
				genMapRed="No_Map_Option"

			fi
			#create checkpoint for the completion of mapping to the host genome

			#Clear the content of the file $mappinglog
			#> $mappinglog

			#tick $checkPoint"/hostMap.chkpt"

			#check that the input file is non-empty before commencing other activities with it
			inFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta"`
                        if [[ $inFileNonEmpty == "yes" ]];
                        then
				#MODULE:NOVEL PREDICTION
				#mirDeep-P prediction
				###########	NOVEL PREDICTION WITH MIRDEEPP	###############################################################################################
				###############################################################################################################################################
				if [[ $mirdeepAnalysis == "yes" ]];
				then
        				#create directory to hold the mirdeep analysis output | or use the same directory and delimit with name : always create directories with the -p flag in case it exists already
					cd $mirDeepP
					logAction "$log" "NOVEL PREDICTION:: Commencing prediction for file $tmp..."
					#refBuilt="yes"$indexDir/$refIndex
					if [[ $refBuilt == "no" ]];
					then
						logAction "$log" "NOVEL PREDICTION:: Reference genome index not yet built. Commencing building..."

						mirReferenceSequenceP=$dataPath"/"$mirReferenceSequence
                                                mirdeepRefName=${mirReferenceSequenceP##*/}
                                                mirdeepRefIndex=${mirdeepRefName%%.*}

						cd $indexDir
                                                prepRefSeq $mirReferenceSequenceP $mirdeepRefIndex
                                                refBuilt="yes"
                                                cd $mirDeepP
					fi

					if [[ $refBuilt == "yes" ]];
                                        then
						logAction "$log" "NOVEL PREDICTION:: Reference genome index available, predicting for $tmp..."
						mirReferenceSequenceP=$dataPath"/"$mirReferenceSequence
						mirdeepRefName=${mirReferenceSequenceP##*/}
                                        	mirdeepRefIndex=${mirdeepRefName%%.*}
						mirOutName=$tmp".mirdeepP"
						mirDeepP $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" $indexDir/$mirdeepRefIndex $mirReferenceSequenceP $dataPath/$mirAnnotationFile $dataPath/$mirChromosomeLength $mirOutName $mirFamilySize $mappinglog $mirDeepP_2 $flankingSequenceLength
					fi

					cd $dataDir/$fname
		        	fi

				#MODULE:PATHOGEN ANALYSIS
				#pathogen analysis
				##############	PATHOGEN REPOSITORY MAPPING HERE	##############################################################################################
				##############################################################################################################################################

				if [[ $microbeAnalysis == "yes" ]];
				#create directory to hold the virus analysis output | or use the same directory and delimit with name
				then
					virMap="nil"
					virMapRed="nil"
                                        subVirMap="nil"
					subVirMapRed="nil"

					microbeReference=$dataPath"/"$microbeReference
					#carry out virus analysis here...first check the index in order to map. If not there, create index
					microbeRefName=${microbeReference##*/}
	        	       		microbeRefIndex=${microbeRefName%%.*}

					checkIndexMicrobeVar=`checkIndex $indexDir $microbeRefIndex`

		               		if [[ $checkIndexMicrobeVar == "no" ]];
	        	       		then
	                	       		logAction "$log" "PATHOGEN MAPPING:: The virus repository reference index does not exist...preparing the virus reference index..."
	                       			cd $indexDir
		                       		prepRefSeq $microbeReference $microbeRefIndex
						cd $dataDir/$fname
	        	               		#cd ../
	               			fi
					#begin mapping to virus here
					logAction "$log" "PATHOGEN MAPPING:: Mapping $tmp to the virus repository..."
	               			BowtieMap $microbeMismatch $indexDir/$microbeRefIndex $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" $dir/$tmp/$tmp".virus.sam" $mappinglog

					if [[ $keepVirMapFiles == "yes" ]];
	                                then
	                                        #move the file to where it will be processed into BAM files
	                                        #create the directory for host map in the maping folder
	                                        virMapFolder=$mappingFilesDir"/Virus"
	                                        mkdir -p $virMapFolder

	                                        cp $dir/$tmp/$tmp".virus.sam" $virMapFolder
	                                fi


					#get only mapped output lines
	        	        	grep -v @ $dir/$tmp/$tmp".virus.sam" > $dir/$tmp/$tmp".virus.txt"
	                		awk '{if($2 == 0 || $2 == 16) {print $0}}' $dir/$tmp/$tmp".virus.txt" > $dir/$tmp/$tmp".virus_mapped.txt"

					#locate virus of interest in the mapping output
					if [ -z "$subjectMicrobe" ]
					then
        					logAction "$log" "PATHOGEN MAPPING:: No subject pathogen selected..."
						subVirMap="-"
                                                subVirMapRed="-"
					else
        					grep $subjectMicrobe $dir/$tmp/$tmp".virus_mapped.txt" > $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt"
					fi

					#grep $subjectMicrobe $dir/$tmp/$tmp".virus_mapped.txt" > $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt"

					virNonEmpty=`fileNonEmpty $dir/$tmp/$tmp".virus_mapped.txt"`
	                        	if [[ $virNonEmpty == "yes" ]];
	                        	then
						#statistics file recording
		                        	virMap=`cut -f 10 $dir/$tmp/$tmp".virus_mapped.txt" | sort | uniq -c | wc -l`
						#virMapRed=`wc -l $dir/$tmp/$tmp".virus_mapped.txt" | awk -F'[^0-9]*' '$0=$2'`
						numUSeq=`cut -f 10 $dir/$tmp/$tmp".virus_mapped.txt" | sort | uniq -c | wc -l`
						numSeq=`cut -f 10 $dir/$tmp/$tmp".virus_mapped.txt" | wc -l`
						numAccession=`cut -f 3 $dir/$tmp/$tmp".virus_mapped.txt" | sort | uniq -c | wc -l`
						#Mock_1_fastq_18-30.seq
						#$dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq"
						#old one: $dir/$tmp/$tmp"_fastq_"$min_length"-"$max_length".seq"
						virMapRed=`awk 'NR==FNR{a[$10]=$0;next} ($1 in a) {print $0}' $dir/$tmp/$tmp".virus_mapped.txt" $dir/$tmp/$tmp"_clip_"$min_length"-"$max_length"_fastq-red.seq" | wc -l`


						cp $dir/$tmp/$tmp".virus_mapped.txt" $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt" $virusOutput

					else
						logAction "$log" "PATHOGEN MAPPING:: No sequence mapped to the virus repository"
						virMap="Empty"
						virMapRed="Empty"
						numUSeq="Empty"
						numSeq="Empty"
						numAccession="Empty"
					fi

					svirNonEmpty=`fileNonEmpty $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt"`
	                        	if [[ $svirNonEmpty == "yes" ]];
	                        	then

						#statistics for subject virus subjectVirus
						subVirMap=`cut -f 10 $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt" | sort | uniq -c | wc -l`

						#subVirMapRed=`cut -f 10 $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt" | wc -l`
						subVirMapRed=`awk 'NR==FNR{a[$10]=$0;next} ($1 in a) {print $0}' $dir/$tmp/$tmp".virus_mapped_"$subjectMicrobe".txt" $dir/$tmp/$tmp"_fastq_"$min_length"-"$max_length".seq" | wc -l | awk -F'[^0-9]*' '$0=$2'`

						#percentSubVirMap=`expr $subVirMap / $nonRedReads`
					else
						logAction "$log" "PATHOGEN MAPPING:: No sequence from $tmp mapped to the virus"
						subVirMap="Empty"
						subVirMapRed="Empty"
					fi
				else
					virMap="No_Map_Option"
					subVirMap="No_Map_Option"
					subVirMapRed="No_Map_Option"
					#$percentGenMap\t$percentSubVirMap
				fi

				#Clear the content of the file $mappinglog
	                        > $mappinglog

				#MODULE:CONSERVED SRNAS
				#Mapping to the miRBase repository
				############	MIRBASE REPOSITORY MAPPING	#################################################################################
				#################################################################################################################################

				if [[ $mirBaseAnalysis == "yes" ]];
	        		#create directory to hold the mirbase analysis output | or use the same directory and delimit with name
		        	then
					mirMap="nil"
	        	        	#carry out mirbase analysis here to identify already known sequence...first check the index in order to map. If not there, create index
					mirBaseReference=$dataPath"/"$mirBaseReference
	                		mirBaseRefName=${mirBaseReference##*/}
		                	mirBaseRefIndex=${mirBaseRefName%%.*}

	        	        	checkIndexMirBaseVar=`checkIndex $indexDir $mirBaseRefIndex`

	                		if [[ $checkIndexMirBaseVar == "no" ]];
		                	then
	        	                	logAction "$log" "MIRBASE MAPPING:: The miRBase repository reference index does not exist...preparing the miRBase reference index..."
	                	        	cd $indexDir
	                        		prepRefSeq $mirBaseReference $mirBaseRefIndex
						cd $dataDir/$fname
		                        	#cd ../
	        	        	fi
	                		#begin mapping to mirBase here
		                	logAction "$log" "MIRBASE MAPPING:: Mapping $tmp to the mirBase repository"
	        	        	BowtieMap $mirBaseMismatch $indexDir/$mirBaseRefIndex $dir/$tmp/$tmp"_clip_allFASTAXX"$min_length"-"$max_length".fasta" $dir/$tmp/$tmp".mirbase.sam" $mappinglog

					if [[ $keepMirMapFiles == "yes" ]];
	                                then
	                                        #move the file to where it will be processed into BAM files
	                                        #create the directory for host map in the maping folder
	                                        mirMapFolder=$mappingFilesDir"/miRBase"
	                                        mkdir -p $mirMapFolder

						#mod: copy only if it is not empty
	                                        miRmapFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp".mirbase.sam"`
        	                                if [[ $miRmapFileNonEmpty == "yes" ]];
                	                        then
                        	                        cp $dir/$tmp/$tmp".mirbase.sam" $mirMapFolder
                                	        fi

	                                	        #cp $dir/$tmp/$tmp".mirbase.sam" $mirMapFolder
	                                fi

					numSubOrgMirMap="nil"
					miRmapFileNonEmpty=`fileNonEmpty $dir/$tmp/$tmp".mirbase.sam"`
					if [[ $miRmapFileNonEmpty == "yes" ]];
                                        then
						#get only mapped output lines
		                		grep -v @ $dir/$tmp/$tmp".mirbase.sam" > $dir/$tmp/$tmp".mirbase.txt"
	        	        		awk '{if($2 == 0 || $2 == 16) {print $0}}' $dir/$tmp/$tmp".mirbase.txt" > $dir/$tmp/$tmp".mirbase_mapped.txt"

						#get family of mapped sequences in mirBase | only if the file is non-zero size
						mirNonEmpty=`fileNonEmpty $dir/$tmp/$tmp".mirbase_mapped.txt"`
	        				if [[ $mirNonEmpty == "yes" ]];
	        				then
							if [[ -z "$subjectOrganismCode" ]];
							then
								numSubOrgMirMap="Option-not-set"
							else
								#numSubOrgMirMap=`grep $subjectOrganismCode $dir/$tmp/$tmp".mirbase_mapped.txt" | cut -f 3 | awk '!seen[$1]++' | grep -v '\-5p' | wc -l`
								numSubOrgMirMap=`grep $subjectOrganismCode $dir/$tmp/$tmp".mirbase_mapped.txt" | cut -f 3 | sort | uniq -c | wc -l`
								cut -f 1,3,10 $dir/$tmp/$tmp".mirbase_mapped.txt" > $dir/$tmp/$tmp"_conserved_mirnas.txt"
							fi
							python $root/"group_mirbase.py" $dir/$tmp/$tmp".mirbase_mapped.txt" $root/"viridiplantae_plants.txt" $dir/$tmp/$tmp".mirbase_mapped_family.txt"
							cp $dir/$tmp/$tmp".mirbase_mapped_family.txt" $mirBaseFamily
							cp $dir/$tmp/$tmp".mirbase_mapped.txt" $mirBaseFamily
							cp $dir/$tmp/$tmp"_conserved_mirnas.txt" $mirBaseFamily
							#mod: here $mirBaseOutput

							#AAGTTCAAGAAAGCTGTGGGA   miR396  Manihot esculenta       Prunus persica  Saccharum officinarum   Populus trichocarpa     Arabidopsis lyrata      Oryza sativa
							#mirbase mapped family output

							#locate miRNA of interest in the mapping output
		                                        grep "${subjectOrganism}" $dir/$tmp/$tmp".mirbase_mapped_family.txt" > $dir/$tmp/$tmp".subject_mapped_mirbase_family.txt"
							#numSubOrgMirMap=`wc -l $dir/$tmp/$tmp".subject_mapped_mirbase_family.txt" | awk -F'[^0-9]*' '$0=$2'`

						else
							numSubOrgMirMap="Empty"
							logAction "$log" "MIRBASE MAPPING:: No sequence from $tmp mapped to mirBase"
						fi

						#statistics file recording
	                        		#mirMap=`cut -f 1 $dir/$tmp/$tmp".mirbase_mapped.txt" | sort | uniq -c | wc -l`
						#mirMap=`cut -f 10 $dir/$tmp/$tmp".mirbase_mapped.txt" | awk -F'[^0-9]*' '$0=$2' | sort | uniq -c | wc -l`
						mirMap=`cut -f 3 $dir/$tmp/$tmp".mirbase_mapped.txt" | awk '{print substr($0,5,length($0))}' | sort | uniq -c | wc -l`
					else
						mirMap="Empty"
						logAction "$log" "MIRBASE MAPPING:: No output from $tmp mapping to mirBase"
					fi
				else
					mirMap="No_Map_Option"
		        	fi

				#Clear the content of the file $mappinglog
                        	> $mappinglog

			fi

			#write the statistics output to the file here
			realSampleName=${FILEMAP[$tmp]}
			statPayLoad="$realSampleName\t$rawReadsLib\t$cleanedReads\t$redReads\t$nonRedReads\t$RfamMap\t$genMapRed\t$genMap\t$virMapRed\t$virMap\t$subVirMapRed\t$subVirMap\t$mirMap\t$numSubOrgMirMap"
			echo -e $statPayLoad >> $readsStat
			echo -e $statPayLoad >> $genStat

			######   create checkpoint file here   #######
			tick $checkPoint"/allRep.chkpt"
			#delete temp input file
			rm $dataDir/"tmpInput.fastq"

			#delete the folder of intermediate files
                	rm -R $dir/$tmp
		done
		logAction "$log" "Ds-Seq:: Completed analysis of $tmp"
		cd $dataDir/
	done
	#tick $checkPoint/$tmp"_hostMap.chkpt"
	tick $checkPoint"/hostMap.chkpt"
	cd $dir/	#move back to the home directory



	#MAKE PLOTS FROM STATISTICS OBTAINED
	####### MAKE PLOTS OF LENGTH DISTRIBUTION       ############################################################
        #Check that the data file exists and non-empty
        dNonEmpty=`fileNonEmpty $smallRNADistro`

        if [[ $dNonEmpty == "yes" ]];
        then
		logAction "$log" "PLOTS:: Making length distribution plot..."
                numLines=`cat $smallRNADistro | wc -l`
                if [[ $numLines -gt 2 ]];
                then
                        if [[ $plotFormat == "pdf" ]];
                        then
                                outFileName=$dir"/"$analysisID"LengthDistributionPlot.pdf"
                        elif [[ $plotFormat == "png" ]];
                        then
                                outFileName=$dir"/"$analysisID"LengthDistributionPlot.png"
                        fi
                        #outFileName=$dir"/LengthDistributionPlot.pdf"
                        Rscript $root/"plotDistro.R" $smallRNADistro $min_length $max_length $numSamples $numRep $outFileName $plotFormat
                        #mv $outFileName $plotsDir    $plotFormat|$plotsDir
                        mv $outFileName $plotsDir
			logAction "$log" "PLOTS:: Length distribution plot successfully created."
                else
                        logAction "$log" "PLOTS:: Length distribution file appears not to have any data to plot"
                fi
        else
                logAction "$log" "PLOTS:: Length distribution file is empty, no plot can be produced"
        fi

	#CALCULATE THE PERCENTAGE READS MAPPING
	dNonEmpty=`fileNonEmpty $readsStat`

        if [[ $dNonEmpty == "yes" ]];
        then
		outCalc="readsMappingPercentages.txt"
		python $root/"sequenceCalculations.py" $readsStat $outCalc

		eFNonEmpty=`fileNonEmpty $outCalc`

        	if [[ $eFNonEmpty == "yes" ]];
        	then
			mv $outCalc $reportsDir
		fi
	else
		logAction "$log" "MAPPING CALCULATION:: Reads statistics file appears empty"
	fi

	#edgeR prediction can only be done after all the replicates have been analyzed
	#the checkpoint file is necessary here
	#create a folder to put the files that will be used for the matrix creation of the edgeR files
	#$checkCheckPoint

	#MODULE:EXPRESSION ANALYSIS
	##########	EXPRESSION ANALYSIS (EDGE-R)	#####################################################################################################
	#####################################################################################################################################################

	checkPoint1=`checkCheckPoint $checkPoint "allRep.chkpt"`
	#modify to check that mapping was done too...before the checkpoint file is checked that it completed
	if [[ $edgeRAnalysis == "yes" ]];
        then
		logAction "$log" "EXPRESSION ANALYSIS:: Commencing EdgeR expression analysis for all samples..."
		#IDEA: What if the user is interested for expression analysis between some samples of interest and not all the samples
		#IDEA: Create a variable that will hold the list of the pairwise sample comparison that can be used instead of all the samples
		#IDEA: To be done on the incremental version

		#get the matrix file first

		if [[ $checkPoint1 == "no" ]];
		then
			logAction "$log" "EXPRESSION ANALYSIS:: The upstream analysis did not complete. This expression analysis stage cannot be initiated. You may need to run the whole analysis again!"
		else
			cd $expDirName
			#get the matrix combination in order to get the sample combination for the calculations, use the array of sample names obtained above for this
			combinationArray=`getArrayCombination "${sampleArray[@]}"`

			if [[ $numRep == 1 ]];
        		then
				#if numRep == 1 then combine the keyfile here
                                logAction "$log" "EXPRESSION ANALYSIS:: No replicates, combining the files..."
                                #key="*_clip_all"$min_length"-"$max_length".txt"
                                key="*_fastq.seq"
                                cat $key | awk '!seen[$0]++' > "all_Key_fastq.seq"
                                #cat $key | cut -d"_" -f 1 | awk '!seen[$0]++' > "all_Key_fastq.seq"
			fi

			for i in ${combinationArray[@]}
			do
        			IFS=_ read sample1 sample2 <<< $i
				#echo "The samples are $sample1 and $sample2"
				#the matrix file calculation for each one is obtained from here
				logAction "$log" "EXPRESSION ANALYSIS:: Generating matrix for combination $i using samples: $sample1 and $sample2 .."

				if [[ $numRep == 1 ]];
                        	then
					logAction "$log" "EXPRESSION ANALYSIS:: Generating expression matrix for no replicates samples..."
					#generateMatrix $sample1 $sample2 $numRep "all_Key_fastq.seq"
					awk 'fname != FILENAME { fname = FILENAME; idx++ } idx == 1 {key[$0] = $0 } idx == 2 {if($1 == key[$1]){ f1[$1] += 1 }} idx == 3 {if($1 == key[$1]){ f2[$1] += 1 }} END {for(seq in key) print seq "\t" f1[seq] "\t" f2[seq] }' "all_Key_fastq.seq" $sample1"_1_clip_"$min_length"-"$max_length"_fastq.seq" $sample2"_1_clip_"$min_length"-"$max_length"_fastq.seq" > $sample1"_"$sample2"_expMatrix.txt"

					#Mock_1_clip_18-25_fastq.seq
                                        #"_clip_"$min_length"-"$max_length"_fastq.seq"
                                        #awk 'NR==FNR{C[$0]=0; next}{for (i=1; i<=NF; i++) if ($i in C) C[$i]++} END{for(i in C) print i "\t" C[i]}' "all_Key_fastq.seq" $sample1"_1_clip_"$min_length"-"$max_length"_fastq.seq" > $sample1".count"
                                        #awk 'NR==FNR{C[$0]=0; next}{for (i=1; i<=NF; i++) if ($i in C) C[$i]++} END{for(i in C) print i "\t" C[i]}' "all_Key_fastq.seq" $sample2"_1_clip_"$min_length"-"$max_length"_fastq.seq" > $sample2".count"

                                        #awk 'NR==FNR{a[$1]=$0;next} $1 in a {print a[$1] "\t" $2}' $sample1".count" $sample2".count" > $sample1"_"$sample2"_expMatrix.txt"
				else
					logAction "$log" "EXPRESSION ANALYSIS:: Generating expression matrix for replicated samples..."
					generateMatrix $sample1 $sample2 $numRep
				fi
			done
			tick $checkPoint"/allMatrix.chkpt"

			#then generate the analysis profile with edgeR ###################################################################
			for i in ${combinationArray[@]}
                        do
				IFS=_ read sample1 sample2 <<< $i
                                logAction "$log" "EXPRESSION ANALYSIS:: Edger calculation for combination $i using samples: $sample1 and $sample2 .."
				#mod: absolute file path???
                                exprOutFile=$sample1"_"$sample2"_expr.txt"
                                matrixFile=$sample1"_"$sample2"_expMatrix.txt"
                                matNonEmpty=`fileNonEmpty $matrixFile`
                                if [[ $matNonEmpty == "yes" ]];
                                then
					#`rExpression $matrixFile $sample1 $sample2 $exprOutFile $numRep1 $numRep2 $pValEdgeR`
                                        #Rscript $edgeRFile $1 $2 $3 $4 $5 $6 $7
					if [[ $numRep == 1 ]];
                                        then
						#do stuff
						logAction "$log" "EXPRESSION ANALYSIS:: Edger calculation (NO REPLICATES) for combination $i using samples: $sample1 and $sample2 .."
						Rscript $root/"edgeR_noRep.r" $matrixFile $sample1 $sample2 $exprOutFile $numRep1 $numRep2 $pValEdgeR
						#Rscript ../pipelineScripts/edgeR_noRep.r M_V_3_matrix.txt "M3" "V3" testtesteMatrixOut.txt 1 1 0.05

						`head -n 1 $sample1"_"$sample2"_expr.txt" > $expDirName/"expHeaderTemp.txt"`
                                                `awk '{print "ID\tSequence\t" $0}' $expDirName/"expHeaderTemp.txt" | sed 's/\"//g' > $expHeader`

                                                #generate the file without the header
                                                `tail -n +2 $sample1"_"$sample2"_expr.txt" > $sample1"_"$sample2"_expr_nH.txt"`

                                                #remove the quotes around the sequence
                                                `sed 's/\"//g' $sample1"_"$sample2"_expr_nH.txt" > $sample1"_"$sample2"_expr_nHQ.txt"`

						#pad the output
                                                awk '{print $0 "\t.\t.\t.\t.\t.\t.\t.\t."}' $sample1"_"$sample2"_expr_nHQ.txt" > $sample1"_"$sample2"_expr_nHQ_noRep_temp.txt"
                                                mv $sample1"_"$sample2"_expr_nHQ_noRep_temp.txt" $sample1"_"$sample2"_expr_nHQ.txt"

					else
						logAction "$log" "EXPRESSION ANALYSIS:: Edger calculation (REPLICATED) for combination $i using samples: $sample1 and $sample2 .."
	                                        Rscript $root/"edgeR.r" $matrixFile $sample1 $sample2 $exprOutFile $numRep1 $numRep2 $pValEdgeR $readCutoff $readSampleCutoff
						#rename the output file due to change in naming convention
						mv $exprOutFile"_DE.txt" $exprOutFile
						mv $exprOutFile"_plotMDS.pdf" $sample1"_"$sample2"_plotMDS.pdf"
						mv $exprOutFile"_plotBCV.pdf" $sample1"_"$sample2"_plotBCV.pdf"

						mv $sample1"_"$sample2"_plotMDS.pdf" $sample1"_"$sample2"_plotBCV.pdf" $plotsDir

	                                        #get the header of the expression file
	                                        `head -n 1 $sample1"_"$sample2"_expr.txt" > $expDirName/"expHeaderTemp.txt"`
	                                        `awk '{print "ID\tSequence\t" $0}' $expDirName/"expHeaderTemp.txt" | sed 's/\"//g' > $expHeader`

	                                        #generate the file without the header
	                                        `tail -n +2 $sample1"_"$sample2"_expr.txt" > $sample1"_"$sample2"_expr_nH.txt"`

	                                        #remove the quotes around the sequence
	                                        `sed 's/\"//g' $sample1"_"$sample2"_expr_nH.txt" > $sample1"_"$sample2"_expr_nHQ.txt"`
					fi

				else
					logAction "$log" "EXPRESSION ANALYSIS:: The matrix file for sample combination ${sample1} and ${sample2} does not exist or is empty, so the expression file cannot be produced"
				fi

			done
			cd $dir/
		fi
	fi
	#$samplesReplicated			---yes or no
	#$numSamples				---number of samples

	#regulation statistics	#############################################################################################
	if [[ $edgeRAnalysis == "yes" ]];
        then
		#create statistics file for the regulation analysis
		#echo -e "Samples\tTotal_output\tDiff_Regulated\tUp-regulated\tDown-regulated" >> $expReport
		checkPointExpr=`checkCheckPoint $checkPoint "allMatrix.chkpt"`
		if [[ $checkPointExpr == "no" ]];
		then
			logAction "$log" "REGULATION STATISTICS:: The expression analysis module did not complete"
		else
			echo -e "\n\nREGULATION ANALYSIS\n" >> $genStat
			echo -e "Samples\tTotal_output\tDiff_Regulated\tUp-regulated\tDown-regulated" >> $genStat
                	combinationArray=`getArrayCombination "${sampleArray[@]}"`

			#upreg=$logFC
			#downreg=$(echo "0 - $logFC" | bc)

	                for i in ${combinationArray[@]}
        	        do
				IFS=_ read sample1 sample2 <<< $i
				exNonEmpty=`fileNonEmpty $expDirName/$sample1"_"$sample2"_expr_nHQ.txt"`
			        if [[ $exNonEmpty == "yes" ]];
        			then
        				totExpr=`cut -f 1 $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" | sort | uniq -c | wc -l`

					python $root/"diffReg.py" $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt" 1 $logFC $pVal 2 5
					python $root/"diffReg.py" $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $expDirName/$sample1"_"$sample2"_expr_nHQ_upreg.txt" 2 $logFC $pVal 2 5
					python $root/"diffReg.py" $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $expDirName/$sample1"_"$sample2"_expr_nHQ_downreg.txt" 3 $logFC $pVal 2 5

					totReg=`wc -l $expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt"`
					totUpReg=`wc -l $expDirName/$sample1"_"$sample2"_expr_nHQ_upreg.txt"`
					totDownReg=`wc -l $expDirName/$sample1"_"$sample2"_expr_nHQ_downreg.txt"`
					echo -e "$i\t$totExpr\t$totReg\t$totUpReg\t$totDownReg\n" >> $genStat
				fi
			done
			echo -e "\n\n" >> $genStat
		fi
	fi

#######		COMPLETION OF THE EDGER EXPRESSION PROFILE ANALYSIS	###########################################################################################################################


	##################################################################################################################################################
	#############	SEQUENCE SEGREGATION	##########################################################################################################
	##################################################################################################################################################

	#create a folder to put the files of the second layer of segregation
	seqSegregatelevel2=$dir"/seqSeg/seg2"
	mkdir -p $seqSegregatelevel2

	#to get the sequence segregation lists

	hostMap=`checkCheckPoint $checkPoint "hostMap.chkpt"`
	if [[ $hostMap == "no" ]];
        then
		logAction "$log" "SEQUENCE SEGREGATION:: The host mapping analysis did not complete. This stage cannot be initiated. You may need to run the host mapping analysis again!"
	else
		cd $seqSegregate
		#use the matrix of samples to combine them
		uniq=($(printf '%s\n' "${sampleArray[@]}" | sort -u))

	        for i in "${uniq[@]}";
        	do
			key=$i"_*"
                	`cat $key | awk '!seen[$1]++' > $i"_all_seq.txt"`
			#mod: copy if it is not empty
			matFileNonEmpty=`fileNonEmpty $i"_all_seq.txt"`
                        if [[ $matFileNonEmpty == "yes" ]];
                        then
				cp $i"_all_seq.txt" $seqSegregatelevel2
			fi
	        done
	fi
	#place a checkpoint for the completion of this phase...it is to be used in the next
	tick $checkPoint"/seg1.chkpt"

	#To obtain all the non-redundant sequences from all samples (to be used when getting the unique sequences to a genotype)
	seg1=`checkCheckPoint $checkPoint "seg1.chkpt"`
	if [[ $seg1 == "no" ]];
        then
		logAction "$log" "SEQUENCE SEGREGATION:: The sequence segregation analysis did not complete. This stage cannot be initiated!"
	else
		#
		allKey="*_all_seq.txt"
		`cat $allKey | awk '!seen[$1]++' > "all_seq.txt"`
	fi

	#to obtain the sequences that are common to all the genotypes
	if [[ $seg1 == "no" ]];
        then
                logAction "$log" "SEQUENCE SEGREGATION:: The sequence segregation analysis did not complete. This stage cannot be initiated!"
        else
                #obtain the common sequencs for a two-sample analysis first
                uniqB=($(printf '%s\n' "${sampleArray[@]}" | sort -u))
		uniqBLen=${#uniqB[@]}
		if [[ $uniqBLen -eq 2 ]];
		then
			F=${uniqB[0]}
			S=${uniqB[1]}
			`awk 'NR==FNR{a[$1]=$0;next} ($1 in a) {print $0}' $seqSegregate/$F"_all_seq.txt" $seqSegregate/$S"_all_seq.txt" > $seqSegregate/"samples_common.txt"`

		#mod: if the file above is empty, there will be no need to identify for others and a conditional statement will be necessary here. The execution is only contingent on common sequences from here
		#obtain for more than two-samples
		elif [[ $uniqBLen -gt 2 ]];
		then
			F=${uniqB[0]}
                        S=${uniqB[1]}
			logAction "$log" "SEQUENCE SEGREGATION:: For the common sequences of more than 2 samples...the first two are $F and $S"
			ruf=`expr $uniqBLen - 2`	#real roof of the iteration based on exclusion of the first two items

			#combine the first two to serve as the base for the rest
			`awk 'NR==FNR{a[$1]=$0;next} ($1 in a) {print $0}' $seqSegregate/$F"_all_seq.txt" $seqSegregate/$S"_all_seq.txt" > $seqSegregate/"temp_common.txt"`
			cp $seqSegregate/"temp_common.txt" $seqSegregate/"${F}_${S}_temp_common.txt"
			#mv $seqSegregate/"${F}_${S}_temp_common.txt" $permDir

			#IDEA: How do you handle when there is no common sequences between the samples
                	for i in "${!uniqB[@]}";
                	do
				if [[ $i -lt $ruf ]];
				then
					aind=`expr $i + 2`
					elem=${uniqB[$aind]}
					logAction "$log" "SEQUENCE SEGREGATION:: For the common sequences...the others are $elem"
					if [[ -z "$elem" ]];
					then
						:
					else
						`awk 'NR==FNR{a[$1]=$0;next} ($1 in a) {print $0}' $seqSegregate/"temp_common.txt" $seqSegregate/$elem"_all_seq.txt" > $seqSegregate/"temp_temp_common.txt"`
						cp $seqSegregate/"temp_temp_common.txt" $seqSegregate/"${elem}_temp_temp_common.txt"
						#mv $seqSegregate/"${elem}_temp_temp_common.txt" $permDir
						mv $seqSegregate/"temp_temp_common.txt" $seqSegregate/"temp_common.txt"
					fi
				fi
				#mv $seqSegregate/"temp_common.txt" $seqSegregate/"samples_common.txt"
                	done
		fi

        fi
	mv $seqSegregate/"temp_common.txt" $seqSegregate/"samples_common.txt"

	#to further segragate the samples and obtain the squences that are unique to each genotype
	#create a folder to put the second level of segregation files from above and move the samples into them
	#move into the folder created above here for the next steps of the analysis
	#use checkpoint to be sure the previous step completed - done

	if [[ $seg1 == "no" ]];
        then
                logAction "$log" "SEQUENCE SEGREGATION:: The sequence segregation analysis did not complete. This stage cannot be initiated!"
        else
                #cd into the created folder
		cd $seqSegregatelevel2

		uniqA=($(printf '%s\n' "${sampleArray[@]}" | sort -u))

		secondArr=(${uniqA[*]})         #copy to a second array with numeric index

		for el in "${!uniqA[@]}"        #loop through main associative array
		do
        		fval=${uniqA[$el]}
        		mcmd=""                 #initialize accumulator
        		for ind in "${!secondArr[@]}";
        		do
                		val=${secondArr[$ind]}
                		if [[ $val == $fval ]];
                		then
                        		#echo -e "$val\t$fval\tNot to be used"
					:
                		else
                        		#echo -e "$val"				#_all_seq.txt
                        		mcmd="$mcmd ${val}_all_seq.txt"     #accumulate command
                		fi
        		done
        		`cat $mcmd | awk '!seen[$1]++' > $fval"_prime.txt"`             #execute command to create new file

			uniqueNonEmpty=`fileNonEmpty $fval"_prime.txt"`
                	if [[ $uniqueNonEmpty == "yes" ]];
                	then
				`awk 'NR==FNR{a[$1]=$0;next} {if($1 in a){next} else {print $0}}' $fval"_prime.txt" $seqSegregate/"all_seq.txt" > $seqSegregate/$fval"_unique.txt"`
			else
				logAction "$log" "SEQUENCE SEGREGATION:: ${fval}_prime.txt is empty"
			fi
		done

		#create checkpoint for this stage
                tick $checkPoint"/seg2.chkpt"
		cd $dir/
        fi


	#MODULE: UNIQUE IDENTIFIER FOR EACH SEQUENCE UNIFORM ACROSS ALL SAMPLES	####################################################################################
	############################################################################################################################################################
	#combining all the mapping to the host to get the genomic loci of all the sequences
        hostAllMap=`checkCheckPoint $checkPoint "hostMap.chkpt"`

        if [[ $hostAllMap == "no" ]];
        then
                logAction "$log" "The combination of all host mapping loci appears not to be completed. This stage cannot be initiated!"
        else
                #
                cd $hostMapDir
		#Make a directory in the permanent directory
		permHostMap=$permDir"/HostMap"
		mkdir -p $permHostMap

                catKey="*_mapped.txt"
                `cat $catKey | awk '!seen[$1$3$4$10]++' > "all_host_mapped.txt"`
		#what exactly are you trying to achieve with this???....all loci that a sequence maps to

		`cat $catKey | awk '!seen[$10]++' | cut -f 10 > "all_host_mapped_ID.txt"`
		#generate unique IDs across all the samples for this analysis
		awk '{print "'"$analysisID"'" "-" NR "\t" $0}' "all_host_mapped_ID.txt" > "all_seq_ID.txt"
		#awk '{print "'"$analysisID"'" "-" NR "\t" $0}' "all_host_mapped_ID.txt" | `cut -f 1,11` > "all_seq_ID.txt"

		#cp "all_seq_ID.txt" $permDir

		#obtain the information on the mapping of each sample with the analysis ID and the genomic loci
		#this will then be grouped into families of genomic loci [sort the chromosome loci...Chromosome10_234898, Chromosome02_874839 etc...]

		uniqA=($(printf '%s\n' "${sampleArray[@]}" | sort -u))

                for i in "${uniqA[@]}";
                do
			OSeqDir=$expOutput"/${i}"
                	mkdir -p $OSeqDir

                        hkey=$i"_*"
			`cat $hkey | awk '!seen[$3$4$10]++' > $i"_all_mapped.txt"`
			#***
			awk 'NR==FNR{a[$10]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $i"_all_mapped.txt" "all_seq_ID.txt" | cut -f 1,3-16 > $i"_all_mapped_ID.txt"
			#TST-1   A-139_25_x3     16      Chromosome07    4641027 255     25M     *       0       0       CCATTAACCGCTCGGGCATCGACCC       IIIIIIIIIIIIIIIIIIIIIIIII       XA:i:0  MD:Z:25 NM:i:0

			python $root/"group_host_locus.py" $i"_all_mapped.txt" $i"_all_mapped_locus.txt"
			#CCATTAACCGCTCGGGCATCGACCC       Chromosome07_4641027_- ,

			#cp $i"_all_mapped.txt" $i"_all_mapped_ID.txt" $i"_all_mapped_locus.txt" $permHostMap
			cp $i"_all_mapped.txt" $permHostMap
			cp $i"_all_mapped_ID.txt" $permHostMap
			cp $i"_all_mapped_locus.txt" $permHostMap

                	firstNonEmpty=`fileNonEmpty $i"_all_mapped_ID.txt"`
                	if [[ $firstNonEmpty == "yes" ]];
                	then
                        	scombinationArray=`getArrayCombination "${sampleArray[@]}"`

                        	for a in ${scombinationArray[@]}
                        	do
                                	IFS=_ read sample1 sample2 <<< $a

                                	if [[ ( $i == $sample1 || $i == $sample2 ) ]];
                                	then

						#correlate expression information with analysisID and sequence
						#***
						awk 'NR==FNR{a[$1]=$0;next} ($11 in a){print $0 "\t" a[$11]}' $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $i"_all_mapped_ID.txt" | cut -f 1,16-28 > $i"_all_mapped_ID_expr_${sample1}${sample2}.txt"
						#TST-3   CTAACAGACCGGTAGACTTGAAC -4.66193157535958       11.7318419614414        0.215790396862234       0.785794988294942       3.23323091519675        0       5       0       6       0       0	0

						#generate data for volcano plot and make the plot
						#cut -f 1,2,3 testHeader.test | awk 'BEGIN{print "header"}1' > testHeaderH.test
						cut -f 1,3,5 $i"_all_mapped_ID_expr_${sample1}${sample2}.txt" | awk 'BEGIN{print "seq\tlogFC\tpVal"}1' > $sample1"_"$sample2"_volcano.txt"
						#if file is not empty, make the plot
						#$sample1"_"$sample2"_volcano.txt"
						vaNonEmpty=`fileNonEmpty $sample1"_"$sample2"_volcano.txt"`
                        			if [[ $vaNonEmpty == "yes" ]];
                        			then
							#Rscript ~/volcanoPlot.R ~/data4Volcano-100.txt testVolcanoScript 1.5 0.05 "pdf"
							volcanoOutput=$sample1"_"$sample2"_volcano_all"
							Rscript $root/"volcanoPlot.R" $sample1"_"$sample2"_volcano.txt" $volcanoOutput $logFC $pVal $plotFormat

							va1NonEmpty=`fileNonEmpty $volcanoOutput"."$plotFormat`
							if [[ $va1NonEmpty == "yes" ]];
                                                	then
								mv $volcanoOutput"."$plotFormat $plotsDir
							fi
						else
							logAction "$log" "VOLCANO PLOT:: Sorry, there seems to be no file to plot for $sample1 and $sample2"
						fi



						#correlate expression information with locus
						awk 'NR==FNR{a[$2]=$0;next} ($1 in a){print a[$1] "\t" $0}' $i"_all_mapped_ID_expr_${sample1}${sample2}.txt" $i"_all_mapped_locus.txt" > $i"_all_mapped_ID_locus_expr_${sample1}${sample2}.txt"
						#CCATTAACCGCTCGGGCATCGACCC       Chromosome07_4641027_- ,

						awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0}' $expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt" $i"_all_mapped_ID_locus_expr_"${sample1}${sample2}".txt" > $i"_all_mapped_ID_locus_expr_"${sample1}${sample2}"_diff.txt"
						#get ones that are differentially expressed too
						#upreg=$logFC
			                        #downreg=`expr 0 - $logFC`
						#downreg=$(echo "0 - $logFC" | bc)
                        			#`awk '{ if ((($3 >= "'"$upreg"'") || ($3 <= "'"$downreg"'")) && ($5 < "'"$pVal"'")) {print $0}}' $i"_all_mapped_ID_locus_expr_${sample1}${sample2}.txt" > $i"_all_mapped_ID_locus_expr_${sample1}${sample2}_diff.txt"`

						mv $i"_all_mapped_ID_locus_expr_${sample1}${sample2}.txt" $i"_all_mapped_expression_${sample1}${sample2}.txt"
						mv $i"_all_mapped_ID_locus_expr_${sample1}${sample2}_diff.txt" $i"_all_mapped_expression_${sample1}${sample2}_diff.txt"
						#mv $i"_all_mapped_expression_${sample1}${sample2}.txt" $i"_all_mapped_expression_${sample1}${sample2}_diff.txt" $OSeqDir
						cp $i"_all_mapped_expression_${sample1}${sample2}.txt" $OSeqDir
						cp $i"_all_mapped_expression_${sample1}${sample2}_diff.txt" $OSeqDir
					#else
					fi
				done
			#else
			fi
		done

		#all genomic loci for all samples
		locKey="*_all_mapped_locus.txt"
		`cat $locKey | awk '!seen[$1]++' > "all_mapped_locus.txt"`
                cd $dir/
        fi
#################################################################################################################################
##########      EXPRESSION PROFILING OF THE SEGREGATED READS AND OTHERS         #################################################
#################################################################################################################################

        #check the existence of the file before attempting to use it
	#common to samples | unique to either samples |
	commonNonEmpty=`fileNonEmpty $seqSegregate/"samples_common.txt"`

        if [[ $commonNonEmpty == "yes" ]];
        then			#loop through the array combination to get the common sequences between the
		commSeqDir=$expOutput"/commonSequences"
		mkdir -p $commSeqDir

		#get the loci of the sequences before looping
		#***
		#to get the common samples analysis ID
                awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0}' $seqSegregate/"samples_common.txt" $hostMapDir/"all_seq_ID.txt" > $seqSegregate/"samp_common_ID.txt"
                awk 'NR==FNR{a[$2]=$0;next} ($1 in a){print a[$1] "\t" $0}' $seqSegregate/"samp_common_ID.txt" $hostMapDir/"all_mapped_locus.txt" | cut -f 1,2,4- > $seqSegregate/"samp_comm_ID_loc.txt"
		#numCommReg=`cut -f 1 $seqSegregate/"samp_common_expr_ID_loc.txt" | sort |  uniq -c | wc -l`
		#begin looping here
		combinationArray=`getArrayCombination "${sampleArray[@]}"`

		echo -e "\n\nREADS SEGREGATION\n" >> $genStat
        	echo -e "Sample\tNum_seqs\tNum_seqs_expression\tNum_seqs_regulated" >> $genStat

                for a in ${combinationArray[@]}
                do
                	IFS=_ read sample1 sample2 <<< $a
			logAction "$log" "UNIQUE ID:: Expression profiling the combination $a for samples $sample1 and $sample2...."
			numComm=`cut -f 1 $seqSegregate/"samples_common.txt" | sort | uniq -c | wc -l`
			awk 'NR==FNR{a[$1]=$0;next} ($1 in a){print $0}' $seqSegregate/"samples_common.txt" $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" > $seqSegregate/"samples_common_expression_${sample1}${sample2}.txt"
			numCommExpr=`cut -f 1 $seqSegregate/"samples_common_expression_${sample1}${sample2}.txt" | sort | uniq -c | wc -l`

			#awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $1 "\t" a[$2]}' $seqSegregate/"samples_common.txt" $hostMapDir/"all_seq_ID.txt" > $seqSegregate/"samp_common_ID.txt"
			#awk 'NR==FNR{a[$2]=$0;next} ($1 in a){print a[$1] "\t" $2}' $seqSegregate/"samp_common_ID.txt" $hostMapDir/"all_mapped_locus.txt" > $seqSegregate/"samp_comm_ID_loc.txt"

			#***
			awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $seqSegregate/"samples_common_expression_${sample1}${sample2}.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,3- > $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID.txt"
			#get the genomic loci for the sequences
			awk 'NR==FNR{a[$2]=$0;next} ($1 in a){print a[$1] "\t" $0}' $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID.txt" $hostMapDir/"all_mapped_locus.txt" > $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_loc.txt"
			#Filter the regulated ones between the two samples 1. ID,  2. Seq,
			upreg=$logFC
			#downreg=`expr 0 - $logFC`
			#downreg=$(echo "0 - $logFC" | bc)
			awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0}' $expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt" $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_loc.txt" > $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_diff.txt"

			#statistics info here
                        #nummber common to both samples
			numCommReg=`cut -f 1 $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_diff.txt" | sort |  uniq -c | wc -l`

			`cat $expHeader $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_loc.txt" > $seqSegregate/"common_exp_${sample1}${sample2}_ID.txt"`
			`cat $expHeader $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_diff.txt" > $seqSegregate/"common_exp_${sample1}${sample2}_ID_diff.txt"`

			#rename the files before moving to the permanent directory
			mv $seqSegregate/"common_exp_${sample1}${sample2}_ID.txt" $seqSegregate/"common_sequences_expression_${sample1}${sample2}.txt"
			mv $seqSegregate/"common_exp_${sample1}${sample2}_ID_diff.txt" $seqSegregate/"common_sequences_expression_${sample1}${sample2}_diff.txt"
			#mv $seqSegregate/"samp_comm_ID_loc.txt" $seqSegregate/"samples_common.txt"

			#$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput
			cp $seqSegregate/"common_sequences_expression_${sample1}${sample2}.txt" $seqSegregate/"common_sequences_expression_${sample1}${sample2}_diff.txt" $commSeqDir
			#cp $seqSegregate/"samples_common.txt" $permDir
			#cp $seqSegregate/"samp_comm.txt" $seqSegregate/"common_expression.txt" $seqSegregate/"common_expression_diff.txt" $permDir

			#statistics info here
			#nummber common to both samples
			#numCommReg=`cut -f 1 $seqSegregate/"samp_common_expr_${sample1}${sample2}_ID_diff.txt" | sort |  uniq -c | wc -l`

			echo -e "$a\t$numComm\t$numCommExpr\t$numCommReg" >> $genStat
		done
	else
		logAction "$log" "UNIQUE ID:: Sorry, there seems to be no sequences common to the two samples"
		numComm=0
		numCommExpr=0
		numCommReg=0
	fi

	#expression profiling of reads that map to mirbase
	#$dir/$tmp/$tmp".mirbase_mapped.txt"   |   #with family $dir/$tmp/$tmp".mirbase_mapped_family.txt"
	#put all mirbase_family in a directory to be used here...they need to be combined per sample as done for the earlier ones

	if [[ $mirBaseAnalysis == "yes" ]];
	then
		cd $mirBaseFamily
		uniqA=($(printf '%s\n' "${sampleArray[@]}" | sort -u))

	        for i in "${uniqA[@]}";
        	do
			#USeqDir=$expOutput"/${i}"
                	#mkdir -p $USeqDir

			mirKey=$i"_*.mirbase_mapped.txt"
			`cat $mirKey | awk '!seen[$10]++' > $i"_all.mirbase_mapped.txt"`
			#$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput
			mirMapFileNonEmpty=`fileNonEmpty $i"_all.mirbase_mapped.txt"`

                        if [[ $mirMapFileNonEmpty == "yes" ]];
                        then
				cp $i"_all.mirbase_mapped.txt" $mirBaseOutput		#this is located in the permanent directory
			#else
			fi

			#`cut -f 10 $i"_all.mirbase_mapped.txt" | sort | uniq -c | wc -l > $i"_mir_useq_mapped.txt"`

			mirFamKey=$i"_*.mirbase_mapped_family.txt"
			`cat $mirFamKey | awk '!seen[$1]++' > $i"_all.mirbase_mapped_family.txt"`

			#cp $i"_all.mirbase_mapped_family.txt" $mirBaseOutput
			#all these actions with expression data presupposes that expressionAnalysis is "yes"
			#make plans for presenting some output too even if the expressionAnalysis is "no"
			miRfamFileNonEmpty=`fileNonEmpty $i"_all.mirbase_mapped_family.txt"`
			if [[ $miRfamFileNonEmpty == "yes" ]];
                        then
				if [[ $edgeRAnalysis == "yes" ]];
		        	then
					mSeqDir=$expOutput"/miRBase/${i}"
	                        	mkdir -p $mSeqDir

					scombinationArray=`getArrayCombination "${sampleArray[@]}"`
	                        	for a in ${scombinationArray[@]}
        	                	do
                	                	IFS=_ read sample1 sample2 <<< $a

                        	        	if [[ ( $i == $sample1 || $i == $sample2 ) ]];
                                		then
							#mirbase family output
							#AAGTTCAAGAAAGCTGTGGGA   miR396  Manihot esculenta       Prunus persica  Saccharum officinarum   Populus trichocarpa     Arabidopsis lyrata      Oryza sativa
							awk 'NR==FNR{a[$1]=$0;next} ($1 in a){print a[$1] "\t" $0}' $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $i"_all.mirbase_mapped_family.txt" > $i"_all.mirbase_mapped_family_expr_${sample1}${sample2}.txt"
							#ID	sequence
							#***++
							awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $i"_all.mirbase_mapped_family_expr_${sample1}${sample2}.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,2,4- > $i"_all.mirbase_mapped_family_expr_${sample1}${sample2}_ID.txt"
							#ensure family is indexed with sequence
							#This step may not be necessary
							#awk 'NR==FNR{a[$2]=$0;next} ($1 in a){print a[$1] "\t" $2 "\t" $3}' $i"_all.mirbase_mapped_family_expr_${sample1}${sample2}_ID.txt" $i"_all.mirbase_mapped_family.txt" > $i"_all.mirbase_mapped_family_expr_ID_${sample1}${sample2}_fam.txt"
							#mv $i"_all.mirbase_mapped_family_expr_ID_${sample1}${sample2}_fam.txt" $i"_mirbase_mapped_fam_expr_${sample1}${sample2}.txt"

							mv $i"_all.mirbase_mapped_family_expr_${sample1}${sample2}_ID.txt" $i"_mirbase_mapped_fam_expr_${sample1}${sample2}.txt"
							#mv $i"_mirbase_mapped_fam_expr_${sample1}${sample2}.txt" $mSeqDir
							cp $i"_mirbase_mapped_fam_expr_${sample1}${sample2}.txt" $mSeqDir
						#else
						fi
					done

				#isolate reads not mapped to miRbase and see if they have expression information
				else 		#move the files withour expression analysis to the destination defined
					#mv ... ...
					#mv $i"_all.mirbase_mapped_family.txt" $mirBaseOutput
					cp $i"_all.mirbase_mapped_family.txt" $mirBaseOutput
				fi
			#else
			fi
			#mv $i"_all.mirbase_mapped.txt" $mirBaseOutput
			cp $i"_all.mirbase_mapped.txt" $mirBaseOutput
		done
		#cd $dir/
	fi

###################################################################################################################################
##############	MIRDEEP-P further analysis: FILTER CONSERVED SEQUENCES AND CORRELATE EXPRESSION PROFILE ###########################
	if [[ $mirdeepAnalysis == "yes" ]];
        then
		hostAllMap2=`checkCheckPoint $checkPoint "hostMap.chkpt"`

		if [[ $hostAllMap2 == "no" ]];
	        then
        	        logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: The combination of all host mapping loci appears not to be completed. This stage cannot be initiated!"
	        else
			cd $mirDeepP_2

			#do something
			uniqA=($(printf '%s\n' "${sampleArray[@]}" | sort -u))

        	        for i in "${uniqA[@]}";
                	do
				#mDeepDir=$expOutput"/NovelPrediction/${i}"
                                #mkdir -p $mDeepDir

				mirKey=$i"_*.mirdeepP.filter_P_prediction"
				`cat $mirKey | awk '!seen[$7]++' > $mirDeepP_2/$i"_U.txt"`
				#Chromosome01    -       A-1585433_21_x71        Chromosome01_1060       32504509..32504529      32504408..32504529      ACCTCCCCTCAAGGGCTTCCG   ACCTCCCCTCAAGGGCTTCCGCGTCGTCATCTCTGGTGGAGTTAGGTAATTGT

				#map the files to the unique analysis ID to be used from here onwards
				#***
				#awk 'NR==FNR{a[$7]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $mirDeepP_2/$i"_U.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,3- > $mirDeepP_2/$i"_U_ID.txt"
				#rename file then move to the permanent directory
				#mv $mirDeepP_2/$i"_U_ID.txt" $mirDeepP_2/$i"_predicted.txt"
				#mv $mirDeepP_2/$i"_predicted.txt" $mirDeepPOutput
				#cp $mirDeepP_2/$i"_predicted.txt" $mirDeepPOutput

				#TST-17346       Chromosome14    +       A-3221536_22_x22        Chromosome14_823        19354767..19354788      19354767..19354857      TAGAGCCAAGAATGACTTGCCG  TAGAGCCAAGAATGACTTGCCGGATAATGATAAGCTA

				f0NonEmpty=`fileNonEmpty $mirDeepP_2/$i"_U.txt"`
				if [[ $f0NonEmpty == "yes" ]];
                                then
					mDeepDir=$expOutput"/NovelPrediction/${i}"
	                                mkdir -p $mDeepDir

					#map the files to the unique analysis ID to be used from here onwards
					#***
	                                awk 'NR==FNR{a[$7]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $mirDeepP_2/$i"_U.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,3- > $mirDeepP_2/$i"_U_ID.txt"
        	                        #rename file then move to the permanent directory
                	                mv $mirDeepP_2/$i"_U_ID.txt" $mirDeepP_2/$i"_predicted.txt"
                        	        #$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput
                                	#mv $mirDeepP_2/$i"_predicted.txt" $mirDeepPOutput
					cp $mirDeepP_2/$i"_predicted.txt" $mirDeepPOutput

                                	#TST-17346       Chromosome14    +       A-3221536_22_x22        Chromosome14_823        19354767..19354788      19354767..19354857      TAGAGCCAAGAATGACTTGCCG  TAG$

					#convert the file to fasta to be used for mapping to miRBase
					`awk '{print ">"$1"\n"$8}' $mirDeepP_2/$i"_predicted.txt" > $mirDeepP_2/$i"_U.fasta"`

					mirBaseRefName=${mirBaseReference##*/}
                	        	mirBaseRefIndex=${mirBaseRefName%%.*}

	                        	checkIndexMirBaseVar=`checkIndex $indexDir $mirBaseRefIndex`

		                        if [[ $checkIndexMirBaseVar == "no" ]];
        		                then
                		                logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: The miRBase repository reference index does not exist...preparing the miRBase reference index..."
                        		        cd $indexDir
                                		prepRefSeq $mirBaseReference $mirBaseRefIndex
		                                cd $mirDeepP_2
        		                fi
                		        #begin mapping to mirBase here
                        		logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: Mapping predicted sequences to the mirBase repository"
	                        	#BowtieMap $mirBaseMismatch $indexDir/$mirBaseRefIndex $mirDeepP_2/$i"_U.fasta" $mirDeepP_2/$i"_mirbase.sam" $mappinglog

					#contine here by using raw bowtie command
					bowtie -a --best -v $mirBaseMismatch $indexDir/$mirBaseRefIndex -f $mirDeepP_2/$i"_U.fasta" --un $mirDeepP_2/$i"_mirbase_unaligned.fasta" -S -p 12 > $mirDeepP_2/$i"_mirbase.sam" 2>> $mappinglog
	        	                #get only mapped output lines
        	        	        grep -v @ $mirDeepP_2/$i"_mirbase.sam" > $mirDeepP_2/$i"_mirbase.txt"
					#check that the file is not empty
					fNonEmpty=`fileNonEmpty $mirDeepP_2/$i"_mirbase.txt"`

					if [[ $fNonEmpty == "yes" ]];
					then
						awk '{if($2 == "0" || $2 == "16") {next} else {print $0}}' $mirDeepP_2/$i"_mirbase.txt" | awk '!seen[$10]++' > $mirDeepP_2/$i"_U_unmapped.txt"
        	                		awk '{if($2 == "0" || $2 == "16") {print $0}}' $mirDeepP_2/$i"_mirbase.txt" > $mirDeepP_2/$i"_mirbase_mapped.txt"
						#awk '{if($2 == "0" || $2 == "16") {print $0}}' $mirDeepP_2/$i"_mirbase.txt" | awk '!seen[$10]++' > $mirDeepP_2/$i"_mirbase_mapped.txt"
						#extract the ones not mapped to mirBase | $i"_U.txt" still in original mirDeepP prediction format while | $i"_U_unmapped.txt" is in mapped file
						awk 'NR==FNR{a[$10]=$0;next} ($7 in a){print $0}' $mirDeepP_2/$i"_U_unmapped.txt" $mirDeepP_2/$i"_U.txt" > $mirDeepP_2/$i"_novel_sRNA.txt"
						cp $mirDeepP_2/$i"_novel_sRNA.txt" $mirDeepPOutput

						#get families of those mapped to mirBase
						#infile plantCode outfile outfile2 outfile3

						mMNonEmpty=`fileNonEmpty $mirDeepP_2/$i"_mirbase_mapped.txt"`
						if [[ $mMNonEmpty == "yes" ]];
                                        	then
							`python $root/"group_mirbase_mirDeepP.py" $mirDeepP_2/$i"_mirbase_mapped.txt" $root/"viridiplantae_plants.txt" $mirDeepP_2/$i"_mirbase_mapped_family.txt" $mirDeepP_2/$i"_mirbase_not-mapped.txt" $mirDeepP_2/$i"_mirbase_statistics.txt"`

							if [[ $edgeRAnalysis == "yes" ]];
						        then
								awk 'NR==FNR{a[$10]=$0;next}{if($7 in a){next} else {print $0}}' $mirDeepP_2/$i"_mirbase_mapped.txt" $mirDeepP_2/$i"_U.txt" > $mirDeepP_2/$i"_U_unmapped.txt"
								#Chromosome01    +       B-930051_22_x198        Chromosome01_104        6405596..6405617        6405549..6405617        TCGGACCAGGCTTCATTCCCCT  GGGAATGTTGTCTGGTTCAAGGTCATATGATGCTTTATAGTTTCTCCTCGGA

								#expression profile of those mapped to mirBase
								awk 'NR==FNR{a[$1]=$0;next} ($10 in a){print a[$10]}' $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $mirDeepP_2/$i"_mirbase_mapped.txt" > $mirDeepP_2/$i"_mirbase_mapped_expression.txt"
								#TCGGACCAGGCTTCATTCCCCT  -1.70230781646542       4.71931703452491        0.165918027259783       0.309015485926636       0.946878060357657       0       400     174     405     198     69      148
								#***+++confirm this line
								awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $mirDeepP_2/$i"_mirbase_mapped_expression.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,3- > $mirDeepP_2/$i"_mirbase_mapped_expression_ID.txt"
								#TST-36287       TTGGACTGAAGGGAGCTCCTT   -2.73700868216621       1.20481715494989        0.0555790440997071      0.130482206644174       1.06279127842295        0       26      34      35      14	3	5

								#get the ones that are regulated
								upreg=$logFC
        	                                                downreg=`expr 0 - $logFC`
	                                                        `awk '{ if ((($3 >= "'"$upreg"'") || ($3 <= "'"$downreg"'")) && ($5 < "'"$pVal"'")) {print $0}}' $i"_mirbase_mapped_expression_ID.txt" > $mirDeepP_2/$i"_mirbase_mapped_expression_ID_diff.txt"`

								#rename the files before moving to permanent directory
		                        			mv $mirDeepP_2/$i"_mirbase_mapped_expression_ID.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_expression.txt"
								mv $mirDeepP_2/$i"_mirbase_mapped_expression_ID_diff.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_expression_diff.txt"
								mv $mirDeepP_2/$i"_mirbase_mapped_family.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_family.txt"
								#$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput | realSampleName=${FILEMAP[$tmp]}
								#mv $mirDeepP_2/$i"_predicted_mirbase_mapped_expression.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_expression_diff.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_family.txt" $mirDeepPOutput
								cp $mirDeepP_2/$i"_predicted_mirbase_mapped_expression.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_expression_diff.txt" $mirDeepP_2/$i"_predicted_mirbase_mapped_family.txt" $mirDeepPOutput
							else
								logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: Expression analysis option not selected"
							fi
						else
							#if there are no files mapped to mirbase, then return the original file
							cp $mirDeepP_2/$i"_U.txt" $mirDeepP_2/$i"_U_unmapped.txt"
							logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: No newly predicted sequence mapped to mirBase"
						fi

						#for the mirBase unmapped
						#expression profile of those not mapped to mirBase
						kMNonEmpty=`fileNonEmpty $mirDeepP_2/$i"_U_unmapped.txt"`
                                                if [[ $kMNonEmpty == "yes" ]];
                                                then
							if [[ $edgeRAnalysis == "yes" ]];
							then
								scombinationArray=`getArrayCombination "${sampleArray[@]}"`

        		                                	for a in ${scombinationArray[@]}
	                	                        	do
                                		                	IFS=_ read sample1 sample2 <<< $a

                                                			if [[ ( $i == $sample1 || $i == $sample2 ) ]];
                                                			then
										#awk 'NR==FNR{a[$1]=$0;next} ($10 in a){print a[$10]}' $expDirName/$sample1"_"$sample2"_expr_nHQ.txt" $mirDeepP_2/$i"_U_unmapped.txt" > $mirDeepP_2/$i"_U_unmapped_expression_${sample1}${sample2}.txt"
										#***
										awk 'NR==FNR{a[$1]=$0;next} ($10 in a){print a[$10]}' $expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt" $mirDeepP_2/$i"_U_unmapped.txt" > $mirDeepP_2/$i"_U_unmapped_expression_${sample1}${sample2}.txt"
										#$expDirName/$sample1"_"$sample2"_expr_nHQ_reg.txt"

										awk 'NR==FNR{a[$1]=$0;next} ($2 in a){print $0 "\t" a[$2]}' $mirDeepP_2/$i"_U_unmapped_expression_${sample1}${sample2}.txt" $hostMapDir/"all_seq_ID.txt" | cut -f 1,3- > $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}.txt"
										#TST-36287       TTGGACTGAAGGGAGCTCCTT   -2.73700868216621       1.20481715494989        0.0555790440997071      0.130482206644174       1.06279127842295        0       26      34      35      14

										#Volcano plot of the newly predicted small RNAs here
										cut -f 1,3,5 $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}.txt" | awk 'BEGIN{print "seq\tlogFC\tpVal"}1' > $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}_volcano.txt"
                                                				#if file is not empty, make the plot
                                                				#$sample1"_"$sample2"_volcano.txt"
                                                				vanNonEmpty=`fileNonEmpty $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}_volcano.txt"`
                                                				if [[ $vanNonEmpty == "yes" ]];
                                                				then
                                                        				#Rscript ~/volcanoPlot.R ~/data4Volcano-100.txt testVolcanoScript 1.5 0.05 "pdf"
                                                        				volcanoOutput=$i"_in_"$sample1"_"$sample2"_volcano"
                                                        				Rscript $root/"volcanoPlot.R" $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}_volcano.txt" $volcanoOutput $logFC $pVal $plotFormat

                                                        				van1NonEmpty=`fileNonEmpty $volcanoOutput"."$plotFormat`
                                                        				if [[ $van1NonEmpty == "yes" ]];
                                                        				then
                                                                				mv $volcanoOutput"."$plotFormat $plotsDir
                                                        				fi
                                                				else
                                                        				logAction "$log" "VOLCANO PLOT:: Sorry, there seems to be no file to plot for $i in $sample1 and $sample2"
                                                				fi

										#Split into up and down-regulated
										#upreg=$logFC
                								#downreg=`expr 0 - $logFC`
                								#`awk '{ if ((($3 >= "'"$upreg"'") || ($3 <= "'"$downreg"'")) && ($5 < "'"$pVal"'")) {print $0}}' $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}.txt" > $mirDeepP_2/$i"_U_unmapped_expression_ID_diff_${sample1}${sample2}.txt"`
										#THESE FILES NEED TO BE MOVED TO THE PERMANENT DIRECTORY
										mv $mirDeepP_2/$i"_U_unmapped.txt" $mirDeepP_2/$i"_predicted_${sample1}${sample2}.txt"
										#mv $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}.txt" $mirDeepP_2/$i"_predicted_unmapped_mirbase_expression_${sample1}${sample2}.txt"
										mv $mirDeepP_2/$i"_U_unmapped_expression_ID_${sample1}${sample2}.txt" $mirDeepP_2/$i"_predicted_expression_${sample1}${sample2}_diff.txt"
										#mv $mirDeepP_2/$i"_U_unmapped_expression_ID_diff_${sample1}${sample2}.txt" $mirDeepP_2/$i"_predicted_unmapped_mirbase_expression_${sample1}${sample2}_diff.txt"
										#$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput
										#mv $mirDeepP_2/$i"_predicted_unmapped_mirbase_expression_${sample1}${sample2}.txt" $mirDeepP_2/$i"_predicted_unmapped_mirbase_expression_${sample1}${sample2}_diff.txt" $mDeepDir
										cp $mirDeepP_2/$i"_predicted_expression_${sample1}${sample2}.txt" $mirDeepP_2/$i"_predicted_${sample1}${sample2}.txt" $mDeepDir
									#else
									fi
								done
							else
								logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: Expression analysis option not selected"
							fi

						else
							logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: No mirbase-unmapped predicted small RNA"
						fi
					else
						logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: File is empty"
					fi
				else
					logAction "$log" "NOVEL (MIRDEEP-P) PREDICTION:: File is empty"
				fi
			cd $dir/
			done
		fi
	fi


	#######		Housecleaning		###################################################################################################################
	#Copy the statistics files into the permanent directory so they dont get deleted when cleaning up
	#$expOutput | $mirDeepPOutput | $mirBaseOutput | $virusOutput

	sRNADistroNonEmpty=`fileNonEmpty $smallRNADistro`
	if [[ $sRNADistroNonEmpty == "yes" ]];
        then
		mv $smallRNADistro $reportsDir
	else
		rm $smallRNADistro
	fi

	fileMapNonEmpty=`fileNonEmpty $fileMap`
        if [[ $fileMapNonEmpty == "yes" ]];
        then
                mv $fileMap $reportsDir
        else
                rm $fileMap
        fi

	readStatNonEmpty=`fileNonEmpty $readsStat`
        if [[ $readStatNonEmpty == "yes" ]];
        then
                mv $readsStat $reportsDir
        else
                rm $readsStat
        fi

	expROutNonEmpty=`fileNonEmpty $expReport`
        if [[ $expROutNonEmpty == "yes" ]];
        then
                mv $expReport $reportsDir
        else
                rm $expReport
        fi

	virOutNonEmpty=`fileNonEmpty $virReport`
        if [[ $virOutNonEmpty == "yes" ]];
        then
                mv $virReport $reportsDir
        else
                rm $virReport
        fi

	genStatOutNonEmpty=`fileNonEmpty $genStat`
        if [[ $genStatOutNonEmpty == "yes" ]];
        then
                mv $genStat $reportsDir
        else
                rm $genStat
        fi


	#convert the host mapping sam files to bam and store in the permanent directory

	#Remove the checkpoint folder
	#Remove the sample folders too
	logAction "$log" "CLEANUP:: Completing analysis, cleaning up intermediate files ..."
	#rm -R $seqSegregate $hostMapDir $expDirName $checkPoint $mirDeepP $mirBaseFamily

	#THIS SHOULD BE LAST...ONLY PRIOR TO DEPLOYMENT,so as to save the time taken to generate index
	#Remove the index directory

	##########################################################################################################################################################

	END_1=$(date +%s)

	DIFF_1=$(( $END_1 - $START_1 ))

	stop=`date +%H:%M:%S`

	echo -e "Script execution completed at $stop == $END_1" >> $log
	echo -e "Script execution completed within $DIFF_1" >> $log

	echo -e "$(($DIFF_1 / 3600)) hours, $((($DIFF_1 / 60) % 60)) minutes and $(($DIFF_1 % 60)) seconds elapsed" >> $log

	echo "Analysis completed!"

	#create file showing analysis completed
	analysisCompleted=$dir"/analysisCompleted.txt"
        rm -f $analysisCompleted
        createLogFile $analysisCompleted
	echo -e "Analysis completed!" >> $analysisCompleted

else
        END_1=$(date +%s)

        DIFF_1=$(( $END_1 - $START_1 ))

        stop=`date +%H:%M:%S`

        echo -e "Script execution completed at $stop == $END_1" >> $log
        echo -e "Script execution completed within $DIFF_1" >> $log
	echo -e "$(($DIFF_1 / 3600)) hours, $((($DIFF_1 / 60) % 60)) minutes and $(($DIFF_1 % 60)) seconds elapsed" >> $log

        echo "Oooops!...Something went wrong"
        echo "Please check the log file for details."

	analysisStopped=$dir"/analysisStopped.txt"
        rm -f $analysisStopped
        createLogFile $analysisStopped
        echo -e "Analysis stopped!\nPlease check the log file for details." >> $analysisStopped

fi
