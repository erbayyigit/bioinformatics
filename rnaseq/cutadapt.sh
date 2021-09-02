cutadapt --info-file=info-file.txt -a AGATCGGAAGAGCACAC -A GATCGTCGGACTGTAGA \
--pair-filter=any --nextseq-trim=20 --pair-filter=both -m=20 \
-o trimmed.R1.fastq -p trimmed.R2.fastq \
r1.fastq r2.fastq > summary_report.txt
