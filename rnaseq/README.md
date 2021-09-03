### Sequences for Adapter Trimmers 
#### R1 adapter to trim NEBNext sRNA or Truseq:  
AGATCGGAAGAGCACAC

#### R2 Adapter to trim NEBNext sRNA:  
GATCGTCGGACTGTAGA    
Longer version:  
GATCGTCGGACTGTAGAACTC  

Feed the above sequences to a trimmer. You can use the test files to test your entries since it contains a read pair with adapter contamination in both R1 and R2.

##  TOOL NOTES
Do not use **"cutadapt"** seem to be corrupting fastq. I got 5 lines records. Ther is chance it might be seqtk but. Needs to be tested.
