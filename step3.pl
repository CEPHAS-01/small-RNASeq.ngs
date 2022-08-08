########usage: perl step3.pl input: control_file experiment_file output_directory
########the program reads input biological replicates (control and experiment)


##usage: perl /biodata/ArabidopsisVirus/software/step3_workflow_3replicates_seq.pl inputfile_control1 inputfile_control2 inputfile_control3 inputfile_exp1 inputfile_exp2 inputfile_exp3 output_file

#!/usr/local/bin/perl
use strict; 


my $file1 = $ARGV[0];   ##control1 file
chomp $file1;
my $file2 = $ARGV[1];	##control2 file
chomp $file2;
my $file3 = $ARGV[2];   ##control3 file
chomp $file3;
my $file4 = $ARGV[3];	##exp1 file
chomp $file4;
my $file5 = $ARGV[4];	##exp2 file
chomp $file5;
my $file6 = $ARGV[5];	##exp3 file
chomp $file6;
my $out_file = $ARGV[6]; #directory of output files
chomp $out_file;

unless ( open(MAC1, $file1))
{
	print"Cannot open file \"$file1\"\n\n";
		exit;
}
my @list1 = <MAC1>;##controllo1
close MAC1;
	
unless ( open(MAC2, $file2))
{
	print"Cannot open file \"$file2\"\n\n";
	exit;
}
my @list2 = <MAC2>;#controllo2
close MAC2;

unless ( open(MAC3, $file3))
{
	print"Cannot open file \"$file3\"\n\n";
	exit;
}
my @list3 = <MAC3>;#controllo3
close MAC3;

unless ( open(MAC4, $file4))
{
	print"Cannot open file \"$file4\"\n\n";
	exit;
}
my @list4 = <MAC4>;#esperimento1
close MAC4;

unless ( open(MAC5, $file5))
{
	print"Cannot open file \"$file5\"\n\n";
	exit;
}
my @list5 = <MAC5>;#esperimento2
close MAC5;

unless ( open(MAC6, $file6))
{
	print"Cannot open file \"$file6\"\n\n";
	exit;
}
my @list6 = <MAC6>;#esperimento3
close MAC6;

my %list1;
my %list2;
my %list3;
my %list4;
my %list5;
my %list6;
my %list;
my %name1;
my %name2;
my %name3;
my %name4;
my %name5;
my %name6;
my @seq1;
my @seq2;
my @seq3;
my @seq4;
my @seq5;
my @seq6;
my $total1;
my $total2;
my $total3;
my $total4;
my $total5;
my $total6;

				
#input file:
#sequenza_lunghezza_rpm-frequenza_log2(rpm-frequenza)_raw-frequenza_nr-seq-totali_name

foreach my $list1 (@list1)   ##controllo1
{	
	chomp $list1;
	my @line = split /\t/, $list1;
	my @line1 = split /_/, $line[0];
	
	
	
		
	my $rawfreq1 = $line1[4];
	
		
	my $key = $line1[0];
	
	$list1{$key}= $rawfreq1;
	
	
	
}

foreach my $list2 (@list2)    ##controllo2
{
	chomp $list2;
	my @line = split /\t/, $list2;
	my @line2 = split /_/, $line[0];
	
	
	my $rawfreq1 = $line2[4];
	
	my $key = $line2[0];

	$list2{$key}= $rawfreq1;	

	
	
}

foreach my $list3 (@list3)    ##controllo3
{
	chomp $list3;
	my @line = split /\t/, $list3;
	my @line3 = split /_/, $line[0];
	
	my $rawfreq1 = $line3[4];
	my $key = $line3[0];
	$list3{$key}= $rawfreq1;

	
	
}
	

foreach my $list4 (@list4)    ##exp1
{
	chomp $list4;
	my @line = split /\t/, $list4;
	my @line4 = split /_/, $line[0];
	
	
	
	my $rawfreq1 = $line4[4];
	
		
	my $key = $line4[0];
	$list4{$key}= $rawfreq1;

	
	
}

foreach my $list5 (@list5)    ##exp2
{
	chomp $list5;
	my @line = split /\t/, $list5;
	my @line5 = split /_/, $line[0];
	
	
	
	my $rawfreq1 = $line5[4];
	
		
	my $key = $line5[0];
	$list5{$key}= $rawfreq1;

	
	
}
foreach my $list6 (@list6)    ##exp3
{
	chomp $list6;
	my @line = split /\t/, $list6;
	my @line6 = split /_/, $line[0];
	
	
	
	my $rawfreq1 = $line6[4];
	
		
	my $key = $line6[0];
	$list6{$key}= $rawfreq1;

	
	
}

############
my $null = 0;
while (my ($key, $value) = each (%list1))   #controllo
{
	
	if ($list2{$key})
	{
		
		
		$list{$key} = $value."_".$list2{$key};			
	}
	else 
	{
		
		$list{$key} = $value."_".$null;			
	}
	my $temp = $list{$key};
	if ($list3{$key})
	{
		
		
		$list{$key} = $temp."_".$list3{$key};			
	}
	else 
	{
		
		$list{$key} = $temp."_".$null;			
	}
	$temp = $list{$key};
	if ($list4{$key})
	{
		
		
		$list{$key} = $temp."_".$list4{$key};			
	}
	else 
	{
		
		$list{$key} = $temp."_".$null;			
	}
	$temp = $list{$key};
	if ($list5{$key})
	{
		
		
		$list{$key} = $temp."_".$list5{$key};			
	}
	else 
	{
		
		$list{$key} = $temp."_".$null;			
	}
	$temp = $list{$key};
	if ($list6{$key})
	{
		
		
		$list{$key} = $temp."_".$list6{$key};			
	}
	else 
	{
		
		$list{$key} = $temp."_".$null;			
	}
}

while (my ($key, $value) = each (%list2))   #controllo1
{
	
	if ($list{$key})
	{
	}
	else
	{
		if ($list1{$key})
		{
			$list{$key} = $list1{$key}."_".$value;	
		}
		else 
		{
		
			$list{$key} = $null."_".$value;			
		}
		my $temp = $list{$key};
		if ($list3{$key})
		{
			$list{$key} = $temp."_".$list3{$key};	
		}
		else 
		{
		
			$list{$key} = $temp."_".$null;			
		}
		$temp = $list{$key};
		if ($list4{$key})
		{
			$list{$key} = $temp."_".$list4{$key};	
		}
		else 
		{
		
			$list{$key} = $temp."_".$null;			
		}
		$temp = $list{$key};
		if ($list5{$key})
		{


			$list{$key} = $temp."_".$list5{$key};			
		}
		else 
		{

			$list{$key} = $temp."_".$null;			
		}
		$temp = $list{$key};
		if ($list6{$key})
		{


			$list{$key} = $temp."_".$list6{$key};			
		}
		else 
		{

			$list{$key} = $temp."_".$null;			
		}

	}
}

while (my ($key, $value) = each (%list3))   #exp
{

	if ($list{$key})
	{
	}
	else
	{
		if ($list4{$key})
		{
			$list{$key} = $null."_".$null."_".$value."_".$list4{$key};	
		}
		else 
		{
		
			$list{$key} = $null."_".$null."_".$value."_".$null;			
		}
		my $temp = $list{$key};
		if ($list5{$key})
		{


			$list{$key} = $temp."_".$list5{$key};			
		}
		else 
		{

			$list{$key} = $temp."_".$null;			
		}
		$temp = $list{$key};
		if ($list6{$key})
		{


			$list{$key} = $temp."_".$list6{$key};			
		}
		else 
		{

			$list{$key} = $temp."_".$null;			
		}		
		
		
	}
}
while (my ($key, $value) = each (%list4))   #exp
{

	if ($list{$key})
	{
	}
	else
	{
		if ($list5{$key})
		{
			$list{$key} = $null."_".$null."_".$null."_".$value."_".$list5{$key};	
		}
		else 
		{
		
			$list{$key} = $null."_".$null."_".$null."_".$value."_".$null;			
		}
		my $temp = $list{$key};
		if ($list6{$key})
		{


			$list{$key} = $temp."_".$list6{$key};			
		}
		else 
		{

			$list{$key} = $temp."_".$null;			
		}

		
		
	}
}	








unless ( open(OUT, ">$out_file"))
{
       print"Cannot open file \"$out_file\"\n\n";
       exit;
}

while (my ($key, $value) = each (%list))



{
	
	my @value = split /_/, $value;

	$value =~s/_/\t/g;

  		print OUT "$key\t$value\n";
		
 		

	
	
}	

close OUT;
exit;