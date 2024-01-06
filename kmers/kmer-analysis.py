from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

''' 
Erbay Yigit. Jan 5, 2023
analyses kmers from a FASTA file
'''
input_f  = 'mstad.fasta'
k = 6
kmerfreq = {}
with open(input_f, 'r') as input_handle:
    for seq_record in SeqIO.parse(input_handle, "fasta"):
        sequence = str(seq_record.seq)
        #this generates overlapping kmers(sliding)  
        for i in range(0, len(sequence)-k+1):
            kmer = (sequence[i:i+k])
            if kmer not in kmerfreq:
                kmerfreq[kmer] = 0
            kmerfreq[kmer] += 1    


# Dictionary as argument to Counter
counter = Counter(kmerfreq)
_keys  = []
_values = []
for key, value in counter.most_common(20):
    print(key, value)
    _keys.append(key)
    _values.append(value)

plt.scatter(_keys, _values)
plt.xticks(rotation='vertical', family='courier')
plt.show()
