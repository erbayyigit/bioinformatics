#!/usr/bin/perl
use 5.30.2;
use warnings;
#EYigit 2021.08.04
# Script: 
# find uniq number of genes obrained from attibutes, after choosing BestREfSeq and mRNA from columns.



my $filename = "GRCh38_latest_genomic_20210804.gff"; 
open(INFILE, "<", $filename) || die "$0: can't open $filename for reading: $!";
open(O_TMP_FH, ">", "BestRefSeq_mRNA.gff") || die "$0: can't open $filename for reading: $!";

my %myGenes;
while(<INFILE>){
	next if /^#/;
	chomp;
	my $line = $_;
	my ($source, $feature, $start, $stop, $attribute)=(split /\t/, $line)[1,2,3,4,8];	
	
	if(($source eq "BestRefSeq") and ($feature eq "mRNA")){
		if ( $attribute =~ m/;gene=(.+?);/){#capture gene name
			$myGenes{$1}=1; #assign gene names to keys
		}		
	}
}

my $size = keys %myGenes;
say $size;
close INFILE;
close O_TMP_FH;



