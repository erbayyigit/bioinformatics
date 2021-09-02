#!/usr/bin/perl
use 5.20.0;
use warnings;
use File::Basename;
use Cwd 'abs_path';

my @files = <*.fastq.gz>; #globe gets list of files in the diredtory

my %dataset_names;
my $r1_file;
my $r2_file;
foreach my $file (@files){
	my $abs_path = abs_path($file);
	my $base = basename($abs_path);
	if ($base =~ m/(.+).(?:1|2).fastq.gz/){. #assumes working with .gz file not .fastq.
		my $datasetName = $1;
		$dataset_names{$datasetName}=1;
	}	
}

for my $k (keys %dataset_names){
	$r1_file = "$k".".1.fastq.gz";
	$r2_file = "$k".".2.fastq.gz";
	system("cutadapt --info-file=$k.'info-file.txt' -a AGATCGGAAGAGCACAC -A GATCGTCGGACTGTAGA --pair-filter=any --nextseq-trim=20 --pair-filter=both -m=20 -o 'clean_'$r1_file -p 'clean_'$r2_file $r1_file $r2_file > 'summary_report_'$k"
	)
}



