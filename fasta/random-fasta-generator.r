#library(seqinr)
Num_Sequences = 100
#Num_Sequences = runif(1, min=2, max=10)
#set.seed(1000) #This can be commented out for less definable patterns

DNA = c ("A", "T", "C", "G")
DNA_probabilities = c(0.2, 0.2, 0.3, 0.3)

for (i in 1:Num_Sequences){
  seq_length = sample(60:300, 1)
  line_header = paste("Chowder_", i, sep="")
  DNA_sequence = paste(sample(DNA, seq_length, 
                              replace = TRUE, 
                              prob = DNA_probabilities),
                              collapse = "")

write.fasta(DNA_sequence, line_header, 
            file.out=paste(Sys.Date(), "_random_", Num_Sequences, ".fasta", sep=""), 
            open = "a", 
            nbchar = 60, 
            as.string = TRUE)
}

