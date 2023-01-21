from Bio import SeqIO
import pprint
import textwrap

my_id_file     = open('sample_files2/fasta_id_records.txt','r')
my_fasta_file  = open('sample_files2/fasta_file.fasta','r')
result_file    = open('sample_files2/results.fasta', "w")

my_dictionary = {} # fasta IDs are keys, value can be anything.
for line in my_id_file:
    my_dictionary[line[:-1]] = 'value'

pprint.pprint(my_dictionary)
for seq_record in SeqIO.parse(my_fasta_file, "fasta"):
    if seq_record.id in my_dictionary:
    	# sequence = str(seq_record.seq)
    	# fasta_record = textwrap.fill(sequence, width=60)
    	# print(f">{seq_record.id}\n{fasta_record}")
    	
    	#use this line to write standard fasta record using SeqIO.write
    	SeqIO.write(seq_record, result_file, "fasta") #default wrap is 60 character
    	#SeqIO.write(seq_record, result_file, "fasta-2line")

my_id_file.close()
my_fasta_file.close()
result_file.close()
