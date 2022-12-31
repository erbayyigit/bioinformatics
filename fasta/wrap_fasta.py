'''
Second option:
Use biopython to parse fasta file, and use Python's textwrap nifty module.
'''
from Bio import SeqIO
import textwrap

for seq_record in SeqIO.parse("input.fasta", "fasta"):
	dna = str(seq_record.seq)
	fasta_record=textwrap.fill(dna, width=30)
	print(">",seq_record.id)
	print(fasta_record)
