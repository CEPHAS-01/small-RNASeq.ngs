#!/usr/local/bin/perl -w

###################################################################################################################
# 	
#	HD_remove_four.pl
#
#	For the HD protocol remove 4nt at the beginning and at the end
#	
# 	Andreas Gisel April 2014
#
###################################################################################################################

use strict;

my $file_in = shift;
my $cut = shift;		# nr. nt's to cut at the beginning and end
my $file_out = shift;

unless ( open(IN, $file_in))
{
  	       print"Cannot open file \"$file_in\"\n\n";
  	       exit;
}

unless ( open(OUT, ">$file_out"))
{
	print"Cannot open file \"$file_out\"\n\n";
	exit;
}

my $counter = 0;
my @seq;

while(<IN>)
{
	my $line = $_;
	chomp $line;
	$counter++;
	
	# @HWIEAS210R_0008:6:1:1011:6994#NNTANT/1
	# CAGACCGGTACACTTGAACATCTCGTATGCCGTC
	# +HWIEAS210R_0008:6:1:1011:6994#NNTANT/1
	# PSVUUW[[Z[SUVVXR_^`_]````^^_^_`BBB

	if ($counter == 1)
	{
		push @seq, $line;
	}
	elsif ($counter == 2)
	{
		if (length($line) > 2*$cut)
		{
			my $seq1 = substr($line, $cut, length($line) - 2*$cut);
			push @seq, $seq1;
		}
		else
		{
			push @seq, "short"
		}
	}
	elsif ($counter == 3)
	{
		push @seq, $line;
	}
	elsif ($counter == 4)
	{
		if (length($line) > 2*$cut)
		{
			my $seq1 = substr($line, $cut, length($line) - 2*$cut);
			push @seq, $seq1;
		}
		else
		{
			push @seq, "short"
		}
		
		$counter = 0;
		if($seq[1] ne 'short')
		{
			foreach my $line (@seq)
			{
				print OUT "$line\n";
			}
		}
		@seq = ();
	}

}
