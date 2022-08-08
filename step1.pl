#!/usr/local/bin/perl -w

###################################################################################################################
# 	
#	
#
#	Reads the raw fastq format file and:removes the adaptor or the whole sequnces without adaptors
#	 
#	
#	- input parameters: Filename adaptor_sequences output_name
#
#	
# 	
#
###################################################################################################################

###########considerare solo no_barcode e lasciare lo stesso nome del file
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
	'h',     40,
	'i',     41,);
	
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
	
my $time0 = time;
my $time1 = time;
my $time = $time1-$time0;
print "Start search - $time\n";

# Adaptor_old: CTGTAGGCACCATCAATCGTA
# Adaptor_new: ATCTCGTATGCCGTCTTCTGCTTGT

# Multiplex barcode sequences
# ACTA
# CGGA
# CTCA
# GAGA
# AGTA
# TCAA
# CGGA
# ACTA
# GTAA
# TACA

my $file_in = $ARGV[0];
my $adap = $ARGV[1];
my $dir = $ARGV[2];
my $out = $ARGV[3];
chomp $out;
		
my %exp;			# keeps barcode in key and experiment name in value


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


my $file_out = $dir."/".$out .".fastq";
unless ( open(FASTQ_OUT, ">$file_out"))	# open unique file for adaptor cleaned sequences
{
  	print"Cannot open file \"$file_out\"\n\n";
  	exit;
}
		
$file_out = $dir."/".$out . "_noadaptor.fastq";
unless ( open(FASTQ_OUT2, ">$file_out"))	# open unique file for adaptor cleaned sequences
{
  	print"Cannot open file \"$file_out\"\n\n";
  	exit;
}

unless ( open(FASTQ, $file_in))
{
  	       print"Cannot open file \"$file_in\"\n\n";
  	       exit;
}

my ($header1, $header2, $seq, $score, $avg_score, $counter, $counter2, $counter3,$count_seq, $count1) = (0,0,0,0,0,0,0,0,0,0);
my $exp;
my %nr_seq;
my $count_N = 0;


my $count = 0;
my %stat;

while(<FASTQ>)
{
				# @HWI-EAS269:2:1:8:382#0/1
				# CGGATACGAACCTCCGACCTTGTGGTCTTTAGGCAC
				# +HWI-EAS269:2:1:8:382#0/1
				# a^\]TJVVOV[_Z_BBBBBBBBBBBBBBBBBBBBBB


	$counter++;
	if ($counter == 1)	# extracts the solexa sequence ID
	{
		$header1 = $_;
		chomp $header1;
		$count++;
		$counter2++;
		
		if ($count % 1000000 == 0)
		{
			$time1 = time;
			$time = $time1-$time0;
			print "$count hits - $time\n";
		}
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
		
		if($seq =~ /N/)
		{
			$count_N++; ####counter per N
		}
		else
		{
			if ($adap)
			{
				my @ret = extract_adapt($seq,$score,$adap);
			
				if($ret[1])
				{
					$seq = $ret[0];
					$score = $ret[1];
					my @score = split//, $score;
					my $score_num = 0;
					foreach my $line (@score)
					{
						$score_num = $score_num+$score{$line};
					}
					$avg_score = $score_num/length($seq);

					my $seq_score = $seq.'_'.$score_num.'_'.$score;
					

					if ($avg_score < 12)
					{
						$counter3++;			#counter per tutte le frequenze  con score <12
					}					
					
					else
					{
					
						# print cleaned sequences to corresponding file
						print FASTQ_OUT "$header1\n$seq\n$header2\n$score\n";
						$stat{"exp"}++;
					}	
				}
				else
				{
					# print to file no_adaptor and count
					print FASTQ_OUT2 "$header1\n$seq\n$header2\n$score\n";
					$stat{"exp-no_ad"}++;
				}
			}
		}
		
	}
	
	
}
$file_out = $dir."/log-seq_clean.txt";
unless ( open(LOG, ">$file_out"))	# stat file
{
	print"Cannot open file \"$file_out\"\n\n";
	exit;
}
print LOG "STAT\n";
print LOG "total sequences: $counter2\navg.score < 12: $counter3 (".100/$counter2*$counter3."%)\nSeq N: $count_N\n";

while (my ($key,$value)=each(%stat))
{
	 print LOG "$key: $value\n";
}
close LOG;
close FASTQ_OUT;
close FASTQ_OUT2;

sub extract_adapt
{
	my $seq = shift;
	my $score = shift;
	my $adap = shift;
	
	my $min_adap = substr($adap,0,6);
	
	if ($seq =~ /$min_adap/)
	{
		my $start = length($`);
		my $insert = substr($seq,0,$start);
		my $insert_score = substr($score,0,$start);
		
		return ($insert,$insert_score);
	}
	else
	{
		return 0;
	}
}




