#!/usr/bin/perl
use 5.30.2;
use warnings;
# yigit@neb.com
# calculate total exon lenght sby gen symbol
# you can use this length to dived coverage
# Note if there same genee in multiple chromsomes, script does not consider it (yet)
# sort output file as "sort -k1 > exons_lengths_per_gene.tab"

my %hash;
while(<>){
    chomp(my $line = $_);
    my ($start,$end, $gene)=(split /\t/, $line)[1,2,3];
    $hash{$gene} += ($end-$start);
}

foreach my $k (keys %hash){
    say "$k\t$hash{$k}";
}
__DATA__
chr1	65418	65433	OR4F5
chr1	65519	65573	OR4F5
chr1	69036	71585	OR4F5
chr1	450739	451678	OR4F29
chr1	685715	686654	OR4F16
chr1	923922	924948	SAMD11
chr1	925921	926013	SAMD11
chr1	930154	930336	SAMD11
chr1	931038	931089	SAMD11
chr1	935771	935896	SAMD11
