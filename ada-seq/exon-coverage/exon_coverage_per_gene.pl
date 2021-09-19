#!/usr/bin/perl
use 5.30.2;
use warnings;
#yigit@neb.com Sept 18, 2021
my %geneHash;

while(<>){
    my $line = $_;
    my ($coverage, $gene)=(split /\t+/, $line)[4,9];
    $geneHash{$gene} += $coverage;
    #say $coverage, "\t" ,$gene;
}
foreach my $k (keys %geneHash){
    say "$k $geneHash{$k}";
}

__DATA__
chr1	923939	923940	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	923940	923941	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	923961	923962	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	923963	923964	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	923972	923973	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	924032	924033	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	924037	924038	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	924047	924048	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	924049	924050	.	1	+	chr1	923922	924948	SAMD11	.	+
chr1	924053	924054	.	1	+	chr1	923922	924948	SAMD11	.	+
