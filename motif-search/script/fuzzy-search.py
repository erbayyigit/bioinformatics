''' 
    #https://github.com/taleinat/fuzzysearch
    --This allows fuzzy search, for example it finds a motif with single 
    mismatch.
    --It reads fasta file directly
    Another alternative was recommended but I have not tried yet:
    https://github.com/Martinsos/edlib
    yigit@neb.com May 21, 20023
    USAGE:
    max_deletions=1, "delete" = skip a character in the sub-sequence
    max_insertions=1, maximum # of insertions ("insert" = skip a character in the sequence)

'''
from fuzzysearch import find_near_matches
from Bio import SeqIO

input_file = 'out.fasta'
fasta_records = SeqIO.parse(input_file,'fasta')

for record in fasta_records:
    name, sequence = record.id, str(record.seq)
    subsequence = 'TGGTGT' # distance = 1
    matches = find_near_matches(subsequence, sequence, 
        max_l_dist=1, 
        max_deletions=0,
        max_insertions=0)
    print("-")
    for counter, item in enumerate(matches):
        print (f"{name}\t{counter}, {item}\t{len(record)}")
    

#print(x[0].start, x[0].matched)
#print(sequence[x[0].start:])

