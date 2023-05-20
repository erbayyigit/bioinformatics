''' 
    #https://github.com/taleinat/fuzzysearch
    --This allows fuzzy search, for example find a motif with single 
    mismatch.
    --It read fasta file
    Alternatives recommended but I have not tried yet:
    https://github.com/Martinsos/edlib
    yigit@neb.com May 18, 20023
'''

# Example below.
# Need to write from file version otherwise need to escape line breaks..
from fuzzysearch import find_near_matches
from Bio import SeqIO

mismatch = 1
input_file = 'inputfile.fasta'
fasta_records = SeqIO.parse(input_file,'fasta')

for record in fasta_records:
    name, sequence = record.id, str(record.seq)
    subsequence = 'AAAACAT' # distance = 1
    matches = find_near_matches(subsequence, sequence, max_l_dist=mismatch)
    [print (name,'\t' ,item) for item in matches]


#print(x[0].start, x[0].matched)
#print(sequence[x[0].start:])

