#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $count = 0; 
my $read_lines = 0; 
my $file = $ARGV[0];


my $peptide_id = ""; 
my @uniprot_id = "";
my @sequence = "";
my @phosphosite = "";
my @full = ""; 


my $start_time = localtime();
my $start_run = time();
my $log = "logfile.txt";

###input file normally A08_phosR.txt
###output file e.g. A08_phosR_v2
###in mysql the symbol and the description are added to A08_phosR_v2,
###which are derived from the uniprot ID mapping server. The table A08_phosR_v3 is than genrated and is nput for R.

open(FILE, "$file") || die "scheisse net geklappt!:$!";

while (my $line = <FILE>)
{
	#chomp $line;
	
	$read_lines++;
	
	if ($read_lines % 1001 == 0)
	{ 
		open(LOG2, '>>'.$log) || die "cannot open the logfile!:$!";
		print STDERR "gelesene zeilen: ".$read_lines."\n"; 
		print STDERR 'Benoetigte Zeit: ' . (localtime(time). " - " . $start_time)."\n"; 
		print LOG2 "gelesene zeilen: ".$read_lines."\n"; 
		print LOG2 'Benoetigte Zeit: ' . (localtime(time). " - " . $start_time)."\n"; 
		my $end_run = time (); 
		my $run_time = $end_run - $start_run; 
		print STDERR "Job took $run_time seconds\n";
		print LOG2 "Job took $run_time seconds\n";
		close LOG2;
	}
	
	#if ($read_lines > 1009) { last; }
	
	my @ar = split("\t",$line);

	if ($read_lines == 1) {
		chomp $line;
		print $line."\n"; 
	}

	if ($read_lines > 1)
	{
		$peptide_id = $ar[0];
		while ( $ar[1] =~ /([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})/g )	{ push(@uniprot_id, $1); }
		while ( $ar[2] =~ /\[.{1,2}\]\.(\w+)\.\[.{1,2}\]/g )	{ push(@sequence,$1) };
		#my @multiple_psites = split(";",$ar[3]);
		while ( $ar[3] =~ /\[(.+)\]/g )	{
			my $match = $1;
			if ($match =~ /\]/g ) {
				my @ar2 = split("]", $match);
				$match = $ar2[0];
			}
			while ( $match =~ /(\w{1}\d+)\(/g ) {
				push(@phosphosite,$1);
			}
		}; #\(\d+\)
		#if ( $line =~ /^GN   Name\=(\S+)\;/g )	{ push(@gene_name, $1) };
		shift @ar;
		shift @ar;
		shift @ar;
		shift @ar;
		

		#if ($seq eq "" and @sequence > 1) {next;}
		$count++;

		my $intensities = join("\t",@ar);
		$intensities =~ s/\r//g;
		$intensities =~ s/\n//g;

		#shift @uniprot_id; 
		#shift @sequence; 
		shift @phosphosite;
			
		print $peptide_id."\t".$uniprot_id[1]."\t".$sequence[1]."\t".join("\|",@phosphosite)."\t".$intensities."\n";
		
		@uniprot_id = "";
		@sequence = "";
		@phosphosite = "";
	}
}
close FILE;
print STDERR "die anzahl der geschriebenen zeilen ist: $count\n";
$count = 0;
$read_lines = 0;

print STDERR 'Benoetigte Zeit: ' . (localtime(time). " - " . $start_time)."\n";
my $end_run = time (); 
my $run_time = $end_run - $start_run; 
print STDERR "Job took $run_time seconds\n";



