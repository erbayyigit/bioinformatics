## Workflow-1
#### PrinSEQ.
- Filter minimum Phred score 25
#### SeqPrep Trimmer
- minimum length 20 nt
- AGATCGGAAGAGCACAC
- GATCGTCGGACTGTAGA
#### FastqMcF
- Use it to synchronize reads from from steps\


#### FastQC Run
- Remove eyr354/355 spike-in contaminants

## Workflow -2
Download GFF file with exons \
Convert GFF file to exons.bed file using `gff_to_exons_bed.pl` script. \
Make sure you sort it by bedtools `sort -i file.bed`
Reduce exons using bedtools mergery groupBy. See examples folder. \
`groupBy -i bedtools_merge_byGroup_test_Input.bed -g 1,4 -c 4,2,3 -o count,min,max |awk -v OFS="\t" '{print $1,$4,$5,$2}'`
