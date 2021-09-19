#!/usr/bin/perl
use 5.30.2;
use warnings;
use File::Basename;
#yigit@neb.com
#use Path::Class;
# SCRIPT=convert bedgraph to bed file format with strad information

my @files = glob "../data_bedGraph/*.bedGraph";
foreach my $file (@files){
	my $dir = dirname($file);
	my $filename = basename ($file, '.bedGraph');
	open(INFH, "<", $file) || die "$0: can't open $filename for reading: $!";
	my $ofile = $filename.".bed";
	open(OFH, ">", $ofile) || die "$0: can't open $filename for reading: $!";
	
	while(<INFH>){
		chomp;
		my $line = $_;
		next if /^#/; 
		my($chr, $zStart, $end, $coverage)= (split /\t/, $line)[0,1,2,3];
		if ($coverage =~ s/(\+)(\d+)/$2/g){
			say OFH join("\t", $chr, $zStart, $end, '.', $coverage,  $1);
		}elsif($coverage =~ s/(\-)(\d+)/$2/g){
			say OFH join("\t", $chr, $zStart, $end, '.', $coverage,  $1);
		}	
	}
}
__DATA__
#chr	0-based-pos	pos	strandcoverage
chr1	14703	14704	-1
chr1	14762	14763	-3
chr1	14769	14770	-2
chr1	14771	14772	-1
chr1	14775	14776	-1
chr1	14801	14802	-1
chr1	14809	14810	+2
chr1	15932	15933	+1
chr1	16096	16097	-2
