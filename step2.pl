#!/usr/local/bin/perl -w

###################################################################################################################
# 	
#	step2.pl
#
#	Reads the fastq format and:
#
#	- reads the quality score of each nucleotide to the header of the fasta file
#	- creates a nr sequence fasta file with sequences whose frequency is > than a threshold and without sequences with 'Ns' (name-running number_length(seq)_freq_flag) flag 1 if sequence with score <12
#	- creates a fastq file with all sequences whose frequency is > than a threshold and without sequences with 'Ns'
#	- calculates the frequency of each sequence and adds it to the header
#	- separates the sequences according the size in different files with nr and full sequences whose frequency is > than threshold in fasta e fastq files
#	- creates a statistics file: input_filename_all_min_max.txt (seq_length_rpm-freq_log2(rpm-freq)_raw-freq_nr-seq-tot_name)
#       input parameters: input_filename.fastq output_dir min_length max_length frequency_threshold
# 	
#	
#
###################################################################################################################


use strict;
use Statistics::Lite qw(:all);



my %score;
my %score1 = (
	';',	-5,
	'<',	-4,
	'=',	-3,
	'>',	-2,
	'?',	-1,
	'@',     0,
	'A',     1,
	'B',     2,
	'C',     3,
	'D',     4,
	'E',	 5,
	'F',     6,
	'G',     7,
	'H',     8,
	'I',     9,
	'J',     10,
	'K',     11,
	'L',     12,
	'M',     13,
	'N',     14,
	'O',     15,
	'P',     16,
	'Q',     17,
	'R',     18,
	'S',     19,
	'T',     20,
	'U',     21,
	'V',     22,
	'W',     23,
	'X',     24,
	'Y',     25,
	'Z',     26,
	'[',     27,
	'\\',    28,
	']',     29,
	'^',     30,
	'_',     31,
	'`',     32,
	'a',     33,
	'b',     34,
	'c',     35,
	'd',     36,
	'e',     37,
	'f',     38,
	'g',     39,
	'h',     40,);
	
my %score2 = (
	'!',0,
	'"',1,
	'#',2,
	'$',3,
	'%',4,
	'&',5,
	'\'',6,
	'(',7,
	')',8,
	'*',9,
	'+',10,
	',',11,
	'-',12,
	'.',13,
	'/',14,
	'0',15,
	'1',16,
	'2',17,
	'3',18,
	'4',19,
	'5',20,
	'6',21,
	'7',22,
	'8',23,
	'9',24,
	':',25,
	';',26,
	'<',27,
	'=',28,
	'>',29,
	'?',30,
	'@',31,
	'A',32,
	'B',33,
	'C',34,
	'D',35,
	'E',36,
	'F',37,
	'G',38,
	'H',39,
	'I',40,
	'J',41,
	'K',42,
	'L',43,
	'M',44,
	'N',45,
	'O',46,
	'P',47,
	'Q',48,
	'R',49,
	'S',50,
);
	
	

my $file_in = $ARGV[0];
chomp $file_in;
my $out_dir = $ARGV[1];   #name of directory of output for statistical analysis
chomp $out_dir;
#my $out_dir = "/home/angelica/fasta";
my $minrange = $ARGV[2];   
chomp $minrange;
my $minrange1 = $minrange -1;
my $maxrange = $ARGV[3];
chomp $maxrange;
my $maxrange1 = $maxrange +1; 
my $freq_threshold = $ARGV[4];
chomp $freq_threshold;

my $time0 = time;
my $time1 = time;
my $time = $time1-$time0;
print "Start search - $time\n";

sleep(2);

my @filename= split /\//,$file_in;
my $name= $filename[$#filename];
$name =~s/.fastq//g;

my @name = split /_/,$name;
my $name1 = $name[0];

unless($name)
{
	$name =~ /(.*).fastq/;
}

#### check the scoring system phred or Illumina 1.3

my $grep_count = `grep -c  "f" $file_in`;
print "GREP: $grep_count\n";
sleep (1);
if ($grep_count == 0)
{
	%score = %score2;
}
else
{
	%score = %score1;
}


unless ( open(FASTQ, $file_in))
{
  print"Cannot open file \"$file_in\"\n\n";
  exit;
}

my $file_out2 = $out_dir."/log_".$name.".txt";
unless ( open(LOG, ">$file_out2"))
{
 print"Cannot open file \"$file_out2\"\n\n";
 exit;
}

my $file_out3 = $out_dir."/".$name."_all".$minrange."-".$maxrange.".txt";
unless ( open(ALL, ">$file_out3"))
{
  print"Cannot open file \"$file_out3\"\n\n";
  exit;
}
my $file_out4 = $out_dir."/".$name."_all".$minrange."-".$maxrange.".fasta";
unless ( open(ALLFASTA, ">$file_out4"))
{
	print"Cannot open file \"$file_out4\"\n\n";
	exit;
}

my $file_out4x = $out_dir."/".$name."_allFASTAXX".$minrange."-".$maxrange.".fasta";
unless (open(ALLFASTAXX, ">$file_out4x"))
{
	print"Cannot open file \"$file_out4x\"\n\n";
	exit;
}

my $total_stat = $out_dir."/".$name."_".$minrange."-".$maxrange."_statistics.txt";
unless ( open(TOTAL, ">>$total_stat"))
{
	print"Cannot open file \"$total_stat\"\n\n";
	exit;
}
my $file_out5 = $out_dir."/".$name."sequence_frequencies".".txt";
unless ( open(ALLFREQ, ">$file_out5"))
{
	print"Cannot open file \"$file_out5\"\n\n";
	exit;
}
my $file_out6 = $out_dir."/".$name."_".$minrange."-".$maxrange.".fastq";
unless ( open(FASTQ_ALL, ">$file_out6"))
{
	print"Cannot open file \"$file_out6\"\n\n";
	exit;
}

my ($header1, $header2, $seq, $score, $avg_score, $counter, $counter2, $counter3,$count_seq, $count1, $count_seq_rpm) = (0,0,0,0,0,0,0,0,0,0,0);
my %nr_seq;
#my @count_sig = 0;
#my %count_sig;

my $ofh = select LOG;
$| = 1;
select $ofh;

my @seq;
my @seq_score;
my @sort_seq;
my %freqhash;
while(<FASTQ>)
{
	# @HWI-EAS269:2:1:8:382#0/1
	# CGGATACGAACCTCCGACCTTGTGGTCTTTAGGCAC
	# +HWI-EAS269:2:1:8:382#0/1
	# a^\]TJVVOV[_Z_BBBBBBBBBBBBBBBBBBBBBB
	$counter++;
	if ($counter == 2) # extracts the sequence
	{		
		$seq = $_;		
		chomp $seq;
		if ($freqhash{$seq})
		{
			$freqhash{$seq}++;
		}
		else
		{
			$freqhash{$seq}= 1;
		}
	}
	elsif ($counter == 4) # extracts the score
	{		
		$counter = 0;		
	}				
}
$counter = 0;
close FASTQ;
unless ( open(FASTQ, $file_in))
{
  print"Cannot open file \"$file_in\"\n\n";
  exit;
}
while(<FASTQ>)
{
	
	# @HWI-EAS269:2:1:8:382#0/1
	# CGGATACGAACCTCCGACCTTGTGGTCTTTAGGCAC
	# +HWI-EAS269:2:1:8:382#0/1
	# a^\]TJVVOV[_Z_BBBBBBBBBBBBBBBBBBBBBB
	$count1++;

	if ($count1 % 1000000 == 0)
	{
		$time1 = time;
		$time = $time1-$time0;
		print LOG "$count1 lines - $time\n";
		print "$count1 lines - $time\n";
		


	}

	$counter++;
	if ($counter == 1)	# extracts the solexa sequence ID
	{
		$header1 = $_;
		chomp $header1;
		$counter2++;
		$count1++;
	}
	elsif ($counter == 2) # extracts the sequence
	{
		
		$seq = $_;		
		chomp $seq;

	}
	elsif ($counter == 3)	# extracts the solexa sequence ID
	{
		$header2 = $_;
		chomp $header2;
	}

	elsif ($counter == 4) # extracts the score
	{
		$score = $_;
		chomp $score;
		
		$counter = 0;		
		my @score = split//, $score;
		my $score_num = 0;
		foreach my $line (@score)
		{
			$score_num = $score_num+$score{$line};
		}
		$avg_score = $score_num/length($seq);

		my $seq_score = $seq.'_'.$score_num.'_'.$score;
		push @seq, $seq_score;
		if ($freqhash{$seq} >= $freq_threshold)
		{
			my $file_out = $out_dir."/".$name."_".length($seq).".fastq";
			unless ( open(FASTQ_OUT2, ">>$file_out"))
			{
				print"Cannot open file \"$file_out2\"\n\n";
				exit;
			}
		
			if (length($seq)<$maxrange1 and length($seq)>$minrange1)	
			{	
			

				print FASTQ_ALL $header1."_".length($seq)."\n".$seq."\n".$header2 ."\n".$score."\n";
			}
		
			print FASTQ_OUT2 $header1."_".length($seq)."\n".$seq."\n".$header2 ."\n".$score."\n";
			close  FASTQ_OUT2;
		}
	}
	
				
}
$counter =0;
$seq = "X";
print LOG "number of seqs $#seq\n";
print "number of seqs $#seq\n";
@sort_seq = sort(@seq);
@seq =();
print LOG "sort done\n";
my $count = 0;
my $count_nr = 0;
	
	
for my $i (0..$#sort_seq)
{
	$counter++;

	if ($counter % 100000 == 0)
	{
		$time1 = time;
		$time = $time1-$time0;
		print "$counter hits - $time\n";
		print LOG "$counter seq - $time\n";
	}

			
	if ($sort_seq[$i] =~ /^($seq)_/)
	{
		++$count;

		$nr_seq{$seq} = "$name-$count_nr\t$count\t$avg_score"; 
								
	}
			
	else
	{				
		@seq =();

		$count = 1;
		++$count_nr;
		@seq = split /_/, $sort_seq[$i];
		$seq = $seq[0];
		my $avg_score = $seq[1]/length($seq[0]);
		$nr_seq{$seq} = "$name-$count_nr\t$count\t$avg_score"; 

	}
			
	

					
}

close FASTQ_ALL;
@sort_seq=();
@seq_score=();
$time1 = time;
$time = $time1-$time0;


my %freq;
my %allfreq;
my $k=0; #running number

my (@count, @count_all);
while (my ($key, $value) = each (%nr_seq))
{
	
	$count[length($key)]++;
	$k++;
	my @info = split /\t/, $value;
	my @score = split / /, $info[2];
	$count_all[length($key)]+=$info[1];
	$allfreq{$key}=$info[1];
	#print "$info[1]\n";
	

	if (length($key)<$maxrange1 and length($key)>$minrange1)
	{

		$count_seq_rpm += $info[1];;    #numero sequenze con lunghezza tra 18 e 26
		if ($info[1]>=$freq_threshold)
		{
		my $file_out = $out_dir."/".$name."_".length($key).".fasta";             #crea file fasta per ogni lunghezza di sequenza tra 18 e 26
		unless ( open(FASTA2, ">>$file_out"))
		{
			print"Cannot open file \"$file_out\"\n\n";
			exit;
		}

		#Tayo: my addition starts from here

		my $file_out0 = $out_dir."/fastaXX".$name."_".length($key).".fasta";
		unless (open(FASTAXX, ">>$file_out0"))
		{
		 print"Cannot open file \"$file_out0\"\n\n";
		 exit;
		}



		#print FASTA2 ">".$name1."-".$k."_".length($key)."_".$info[1]."_".$count_sig{$key}."\n$key\n";
						
		#print ALLFASTA ">".$name1."-".$k."_".length($key)."_".$info[1]."_".$count_sig{$key}."\n$key\n";
		print FASTA2">".$name1."-".$k."_".length($key)."_".$info[1]."\n$key\n"; 
		print FASTAXX">".$name1."-".$k."_".length($key)."_x".$info[1]."\n$key\n"; #Tayo: another addition
 
		#if ($info[1]>=$freq_threshold)
		#{
			print ALLFASTA ">".$name1."-".$k."_".length($key)."_".$info[1]."\n$key\n";
			print ALLFASTAXX ">".$name1."-".$k."_".length($key)."_x".$info[1]."\n$key\n";
		#}
		######################################
		#nome-running number_length(seq)_freq es MAC-17-1436_18_245_0 (nome MAC-17 running numb 1436)
		#flag if sequence with score <12			
		}
		$freq{length($key)}{$info[1]}++;	
	}
		
}

my $countfreq=0;
my @sorted = sort {$allfreq{$b} <=> $allfreq{$a}} keys %allfreq;
foreach my $sorted (@sorted)
{
	
	
	if ($allfreq{$sorted}>= $freq_threshold)
	{
		print ALLFREQ "$sorted\t$allfreq{$sorted}\n";

	}
	else
	{
		$countfreq++;  #reads con freq< threshold
	}



}

print LOG "length\tnumber of nr sequences\tnumber of sequences\n";
for (my $i=0;$i<=$#count;$i++)
{
	if($count[$i])
	{
		my $file_out = $out_dir."/".$name."_".$i."_statistics.txt";

		my $file_hist = $out_dir."/".$name."_".$i."_histogram.txt";

		
		
		print LOG "$i\t$count[$i]\t$count_all[$i]\n";

		if ($i<$maxrange1 and $i>$minrange1)
		{	


			unless ( open(FASTA3, ">>$file_out"))
			{
				print"Cannot open file \"$file_out\"\n\n";
				exit;
			}
			unless ( open(HIST, ">>$file_hist"))
			{
				print"Cannot open file \"$file_hist\"\n\n";
				exit;
			}
			
			print FASTA3 "Number of nr sequences size\t$i\t$count[$i]\n";
			print FASTA3 "Number of sequence size\t$i\t$count_all[$i]\n";
			print FASTA3 "Number of sequence size ".$minrange."-".$maxrange."\t$count_seq_rpm\n";
			print TOTAL "$i\t$count[$i]\t$count_all[$i]\n";
			foreach my $key (sort  { $a <=> $b } keys %{$freq{$i}})
			{
				print HIST "$key\t$freq{$i}{$key}\n";
				
			}
			 
		}
	}
}

print LOG "number of reads with freq less than threshold: $countfreq\n";
while (my ($key, $value) = each (%nr_seq))
{
	#if (length($key)<$maxrange1 and length($key)>$minrange1)
	if (length($key)<$maxrange1 and length($key)>$minrange1 and $freqhash{$key}>=$freq_threshold)
	{
		
		my @info = split /\t/, $value;
		my $norm_freq = ($info[1]/$count_seq_rpm)*1000000;
		$norm_freq += 1;	#shift norm_freq by 1 to avoid log(0)
		my $log_norm_freq= log($norm_freq)/log(2);
		print ALL $key."_".length($key)."_".$norm_freq."_".$log_norm_freq."_".$info[1]."_".$count_seq_rpm."_".$info[0]."\n"; #sequenza_lunghezza_rpm-frequenza_log2(rpm-frequenza)_raw-frequenza_nr-seq-totali_name
		
	}
}

$time1 = time;
$time = $time1-$time0;
print LOG "DONE - $time\n";
close FASTA2;
close FASTAXX; #Tayo: addition
close FASTA3;
close TOTAL;
close ALLFASTA;
close ALLFASTAXX; #Tayo: addition
close ALL;
close LOG;
close HIST;
close ALLFREQ;

exit;


