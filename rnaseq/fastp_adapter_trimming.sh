# for small RNA you can use the sequences in the command. It works fing.
# Do not print information file (which is not in the command). It takes too much space.
fastp -i hek_RAND_60m.R1.fastq.gz -o clean_hek_RAND.R1.fastq.gz \
-I hek_RAND_60m.R2.fastq.gz -O clean_hek_RAND.R2.fastq.gz --length_required=25 \
--adapter_sequence=AGATCGGAAGAGCACAC \
--adapter_sequence_r2=GATCGTCGGACTGTAGA 2> reportRAND.txt

