#!/usr/bin/perl
use 5.30.2;
use warnings;
# Erbay Yigit, yigit@neb.com
# Use this file to generage BED regions that you can extract from Zhiyi'start
#bedgraphs to upload teh ucsc browser.
say "#Chr\tStart\tEnd\tSymbol";
while(<>){
	next if /^#/;
	my $line = $_;
	chomp $line;
	my($chr, $strand,  $txStart, $txEnd, $eStart, $eEnd, $gene)=(split /\s+/, $line)[0,1,2,3,5,6,7];
	my $start = ($txStart-1000);
	my $end = ($txEnd+1000);
	my @triple_array = join("\t",$chr, $start, $end, $gene);
	foreach my $item (@triple_array){
		say "$item";
	}
}
__DATA__
#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	name2
chr11	-	1584341	1585283	1	1584341,	1585283,	KRTAP5-1
chr11	-	1597176	1598294	1	1597176,	1598294,	KRTAP5-2
chr11	-	1607564	1608463	1	1607564,	1608463,	KRTAP5-3
chr11	-	1620957	1622138	1	1620957,	1622138,	KRTAP5-4
chr11	+	1629802	1630930	1	1629802,	1630930,	KRTAP5-5
chr11	+	1665598	1667856	3	1665598,1666434,	1665764,1666543,	FAM99A
chr11	-	1683269	1685629	3	1683269,1684193,1685437,	1684034,1684302,1685629,	FAM99B
