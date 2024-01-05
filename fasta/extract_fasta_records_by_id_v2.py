from Bio import SeqIO

ids_file = 'ids.txt'
input_f  = 'input_file.fasta'
output_f = 'result.fasta'

my_dictionary = {}
with open(ids_file, 'r') as f:
    for line in f:
        my_dictionary[line[:-1]]='value'

with open(input_f, 'r') as input_handle, open(output_f, 'w') as output_handle:
    for seq_record in SeqIO.parse(input_handle, "fasta"):
        if seq_record.id in my_dictionary:
            SeqIO.write(seq_record, output_handle, "fasta")

