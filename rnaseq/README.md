# ANALYSIS
#### 1. Downsample
Using ```seqkit sample```
#### 2. Adapter trim
Using ```fastp -i hek_RAND_60m.R1.fastq.gz -o clean_hek_RAND.R1.fastq.gz \
-I hek_RAND_60m.R2.fastq.gz -O clean_hek_RAND.R2.fastq.gz --length_required=20 \
--adapter_sequence=AGATCGGAAGAGCACAC \
--adapter_sequence_r2=GATCGTCGGACTGTAGA 2> reportRAND.txt ``` 
small RNA adapters
#### 3. Leading and Trailing base removal such as random hexamers
Use ```Flexbar``` on Galaxy


### Sequences for Adapter Trimmers 
#### R1 adapter to trim NEBNext sRNA or Truseq:  
```AGATCGGAAGAGCACAC```

#### R2 Adapter to trim NEBNext sRNA:  
```GATCGTCGGACTGTAGA```\
Longer version: ```GATCGTCGGACTGTAGAACTC```

Feed the above sequences to a trimmer. You can use the test files to test your entries since it contains a read pair with adapter contamination in both R1 and R2.

##  TOOL NOTES
Do not use ```cutadapt``` seem to be corrupting fastq. I got 5 lines records. Ther is chance it might be seqtk or sth unrelated. Needs to be tested. Also ```cutadapt``` is not recommneded by experts --too slow and not high fidelity enough.
