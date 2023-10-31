from Bio import SeqIO
import Bio.Data.IUPACData as bdi
from fuzzysearch import find_near_matches
from itertools import product
''' 
    # https://github.com/taleinat/fuzzysearch
    # 
    --It reads fasta file using Biopython
    Another alternative was recommended but I have not tried yet:
    https://github.com/Martinsos/edlib
    yigit@neb.com Oct 30, 20023
    USAGE:
    max_deletions=1, "delete" = skip a character in the sub-sequence
    max_insertions=1, maximum # of insertions ("insert" = skip a character in the sequence)
    
    iupac_symbols = {
    'R' : ['G','A'], 
    'Y' : ['T','C'],
    'M' : ['A','C'],
    'K' : ['G','T'],
    'S' : ['G','C'],
    'W' : ['A','T'],
    'H' : ['A','C','T'],
    'B' : ['G','T','C'],
    'V' : ['G','C','A'],
    'D' : ['G','A','T'],
    'N' : ['G','A','T','C']
}
    USAGE:
    --This works on multirecord FASTA file
    --This can sue iupac symbols in motif sequence, for example H, N, etc. It will iterate
    for each degenerate position.


'''
input_file = 'outputfile.fasta'
motif    = 'AAAAHAT' #motif
max_l_dist_     = 1 #determines distance, 1 is single mismatch
max_deletions_  = 0
max_insertions_ = 0
fasta_records = SeqIO.parse(input_file,'fasta')


def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = bdi.ambiguous_dna_values
   r = []
   for i in product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r
ambiguous_motif = extend_ambiguous_dna(motif)


for record in fasta_records:
    name, sequence = record.id, str(record.seq)    
    for subsequence in ambiguous_motif:
        print(f"Motif {subsequence} in '>{name}':")
        
        matches = find_near_matches(subsequence, sequence, 
            max_l_dist     = max_l_dist_, 
            max_deletions  = max_deletions_, 
            max_insertions = max_insertions_
        ) 
        
        for counter, item in enumerate(matches, 1):
            print (f"\t{counter}\t{item}\t{len(record)}")

    print()



