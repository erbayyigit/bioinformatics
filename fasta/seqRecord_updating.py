from Bio import Seq, SeqIO, SeqRecord
'''
    Following code allows modification of SeqRecord instead of manual printing
    yigit@neb.com. Feb 12, 2023.
'''
input_file   = 'input.fasta'
output_file  = 'output.fasta'

with open(input_file, 'r') as in_handle, open(output_file, 'w') as out_handle:
    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence = str(seq_record.seq)
        trimmed_sequence = sequence[15:-15] #manupulate sequence
        modified_record= SeqRecord.SeqRecord(id=seq_record.id, \
            seq=Seq.Seq(trimmed_sequence), description='description')
        SeqIO.write(modified_record, out_handle, "fasta")
