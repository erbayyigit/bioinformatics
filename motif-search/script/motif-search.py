from Bio import SeqIO
import collections
import re
import os

#this can be very very useful for MazF motif search.

input_file = '/Users/yigit/Downloads/motif-search/data/inputfile.fa'
fasta_sequences = SeqIO.parse(input_file,'fasta')

m_list =[]
num_sequences = 0
for fasta in fasta_sequences:
    num_sequences+=1
    name, sequence = fasta.id, str(fasta.seq)
    matches = re.finditer(r"P[A-Z]{2}P..", sequence)
    #r"P[A-Z]{2}P.." i regex finds like "P..P.." peptides, nothing tricky about it.
    for m in matches:
        match_id = f'{m.start()}_{m.end()}_{m.group()}'
        m_list.append(match_id)
        print(name, m.start(), m.end(), m.group())        
print("total sequences",num_sequences)
print("frequency",collections.Counter(m_list))
print(num_sequences)

