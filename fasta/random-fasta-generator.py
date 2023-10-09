"""
    yigit@neb.com 2023.10.08
    Generates random fasta file.
    Number of records given by user.
    Record length range is giving by user as in k.
    This script performs well as expected.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from random import choices

DNA = ['A', 'T', 'C', 'G']
prob = [0.25, 0.25, 0.25, 0.25] #probability of each base
number_of_records = 100 #number of fasta records

with open("outputfile.fasta", "w") as output_handle:
    for i in range(number_of_records):
        dna_seq = choices(DNA, prob, k=random.randint(100, 500))
        dna_seq = ''.join(dna_seq)
        #Create a SeqRecord object
        record = SeqRecord(Seq(dna_seq),
            id="Chowder_"+str(i+1),
            name="Chow-Chow",
            description="random sequence")
        SeqIO.write(record, output_handle, "fasta")

