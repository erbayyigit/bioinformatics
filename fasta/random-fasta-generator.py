"""
    Generates random fasta file.
    Number of records and length of each record is given by user.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from random import choices


DNA = ('A', 'T', 'C', 'G')
prob = [0.1, 0.1, 0.40, 0.40] #probability of each base
number_of_records = 100 #number of fasta records


with open("outputfile.fasta", "w") as output_handle:
    for i in range(number_of_records):
        dna_seq = choices(DNA, prob, k=random.randint(100, 500))
        dna_str = ''.join(dna_seq)
        record = SeqRecord(Seq(dna_str), #Create a SeqRecord object
            id=f"Chowder_{i+1}",
            name="Chow-Chow",
            description="random sequence")
        SeqIO.write(record, output_handle, "fasta")

