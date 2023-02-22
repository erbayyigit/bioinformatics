from Bio import SeqIO
'''
    yigit@neb.com, 20230221
    converts multiline fasta file to 2line fasta file
'''
infile_name  ='input.fasta'
outfile_name ='output_2line.fasta'

with open(infile_name, 'r') as inFH, open(outfile_name, 'w') as outFH:
    for seq_record in SeqIO.parse(inFH, 'fasta'):
        SeqIO.write(seq_record, outFH, "fasta-2line")

        
        
