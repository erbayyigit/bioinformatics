'''
yigit@neb.com 20220223
'''
from Bio import SeqIO
import re

#following code searches glycosylation motifs in a protein
with open('protein.fasta', 'r') as infile:
    for seq_record in SeqIO.parse(infile, "fasta"):
        #for overlapping ones use look ahead ?=
        #remove ?= for non-overlapping matches
        for m in re.finditer(r"(?=(N[^P][ST][^P]))", str(seq_record.seq)):
            print(m.start(0), m.end(1), m.group(1))

 
