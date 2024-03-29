#!/usr/bin/perl 


use warnings;
use strict;
use Getopt::Std;

######################################### USAGE ################################

my $usage=
"$0 file_SAM_result file_solexa_seq file_chromosome

This is a converter which changes SAM output into Blast format.
The input includes three files: a SAM result file, a fasta file 
consisting of small Reads and a chromosome fasta file. It outputs 
the alignments in blast_parsed format.

It does not matter whether the file_SAM_result file has a header.

file_solexa_seq likes:

>AtFlower100010_x2
AAGGAGATTCTTTCAGTCCAG

file_chromosome contains chromosome seq in fasta format

";


####################################### INPUT FILES ############################

my $file_SAM_result=shift or die $usage;
my $file_short_seq=shift or die $usage;
my $file_chromosome_seq=shift or die $usage;


##################################### GLOBAL VARIBALES #########################

my %short_seq_length=();
my %chromosome_length=();


######################################### MAIN ################################# 

#get the short sequence id and its length
sequence_length($file_short_seq,\%short_seq_length);

#get the chromosome sequence id and its length
sequence_length($file_chromosome_seq,\%chromosome_length);

#convert bowtie result format to blast format;
change_format($file_SAM_result);

exit;


##################################### SUBROUTINES ##############################

sub sequence_length{
    my ($file,$hash) = @_;
    my ($id, $desc, $sequence, $seq_length) = ();

    open (FASTA, "<$file") or die "can not open $$file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $$hash{$id}  = length $sequence;
		    $id         = $1;
		    $desc       = $2;
		    $sequence   = "";
		    next;
                }
		$sequence .= $_;
            }
        }
    }
    $seq_length=length($sequence);
    $$hash{$id} = $seq_length;
    close FASTA;
}





sub change_format{
	#Change Bowtie format into blast format
	my $file=shift @_;
	open(FILE,"<$file")||die"can not open the bowtie result file:$!\n";
  #open(BLASTOUT,">blastout")||die"can not create the blastout file:$!\n";
	
	while(<FILE>){
		chomp;
	#Parse the main part of SAM file	
		if($_=~/^[0-9a-zA-Z]/){
			my @tmp=split("\t",$_);
		#Clean the reads ID
			if($tmp[2] ne '*'){

				my @tmp1=split(" ",$tmp[0]);
				print  "$tmp1[0]"."\t"."$short_seq_length{$tmp1[0]}"."\t"."1".'..'."$short_seq_length{$tmp1[0]}"."\t"."$tmp[2]"."\t"."$chromosome_length{$tmp[2]}"."\t";
				if($tmp[1] == 0){
					my $seq_end=$tmp[3] + $short_seq_length{$tmp1[0]};
					my $seq_bg=$tmp[3] + 1;
					print  "$seq_bg".'..'."$seq_end"."\t"."1e-04"."\t"."1.00"."\t"."42.1"."\t"."Plus / Plus"."\n";
				}
				if($tmp[1] == 16){
					my $seq_end=$chromosome_length{$tmp[2]} - $tmp[3];
					my $seq_bg=$seq_end - $short_seq_length{$tmp1[0]} + 1;
					print  "$seq_bg".'..'."$seq_end"."\t"."1e-04"."\t"."1.00"."\t"."42.1"."\t"."Plus / Minus"."\n";
				}
			}
		}
	}
	
#	close BLASTOUT;

}



