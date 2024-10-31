import argparse
import os
from Bio import SeqIO
import Bio.Data.IUPACData as bdi
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from fuzzysearch import find_near_matches
from itertools import product

'''
    yigit@neb.com
    Nov 5, 2023. Updated on Oct 31, 2024
    Credits: Recoding part of the program was written by Vladimir.
    USAGE:
        ->This script will accept a FASTA file (can be multirecord) as an input. Use DNA bases ie 'T' instead of 'U'
        ->It will recode each FASTA record sequence without changing amino acid sequence
        ->It will search motif in these sequence (will accept 'ambiguous_dna_values" such ash AAAHAT)
        ->Fuzzy search will search the motif in each FASTA record.
        ->Mismatch can be indicated at 'max_l_dist_'. Zero is perfect match.

    REFERENCES:
        https://github.com/taleinat/fuzzysearch
        max_deletions =1, "delete" = skip a character in the sub-sequence
        max_insertions=1, maximum # of insertions ("insert" = skip a character in the sequence)

'''
fasta_file_name     = 'rna_or_ivt_template_references.fasta' #input filename
out_fasta_filename  = f"_recoded_{fasta_file_name}"
outfasta = open(out_fasta_filename, 'w+')       

inpu_file_basename = 'rna_or_ivt_template_references.fasta'.rstrip('.fasta') 
out_match_results  = f"_match_results_{inpu_file_basename}.txt"
outmatch = open(out_match_results , 'w+')       


#fuzzy search parameter
motif           = 'GAAAT'       #'AAAAHAT' #motif
max_l_dist_     = 0 #determines distance, 1 is single mismatch
max_deletions_  = 0
max_insertions_ = 0


codon_table = { 
    #this is Ecoli codon table
    '*': {'TAA': 0.64, 'TAG': 0.07, 'TGA': 0.29},
    'A': {'GCA': 0.21, 'GCC': 0.27, 'GCG': 0.36, 'GCT': 0.16},
    'C': {'TGC': 0.56, 'TGT': 0.44},
    'D': {'GAC': 0.37, 'GAT': 0.63},
    'E': {'GAA': 0.69, 'GAG': 0.31},
    'F': {'TTC': 0.43, 'TTT': 0.57},
    'G': {'GGA': 0.11, 'GGC': 0.41, 'GGG': 0.15, 'GGT': 0.34},
    'H': {'CAC': 0.43, 'CAT': 0.57},
    'I': {'ATA': 0.07, 'ATC': 0.42, 'ATT': 0.51},
    'K': {'AAA': 0.76, 'AAG': 0.24},
    'L': {'CTA': 0.04, 'CTC': 0.1, 'CTG': 0.5, 'CTT': 0.1, 'TTA': 0.13, 'TTG': 0.13},
    'M': {'ATG': 1.0},
    'N': {'AAC': 0.55, 'AAT': 0.45},
    'P': {'CCA': 0.19, 'CCC': 0.12, 'CCG': 0.53, 'CCT': 0.16},
    'Q': {'CAA': 0.35, 'CAG': 0.65},
    'R': {'AGA': 0.04, 'AGG': 0.02, 'CGA': 0.06, 'CGC': 0.4, 'CGG': 0.1, 'CGT': 0.38},
    'S': {'AGC': 0.28, 'AGT': 0.15, 'TCA': 0.12, 'TCC': 0.15, 'TCG': 0.15, 'TCT': 0.15},
    'T': {'ACA': 0.13, 'ACC': 0.44, 'ACG': 0.27, 'ACT': 0.16},
    'V': {'GTA': 0.15, 'GTC': 0.22, 'GTG': 0.37, 'GTT': 0.26},
    'W': {'TGG': 1.0},
    'Y': {'TAC': 0.43, 'TAT': 0.57}
}

wt_seqence = ''
record_list = []

for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
    record_list.append(seq_record) #add wt fasta record
    
    sequence = str(seq_record.seq)
    for base in ['A', 'T', 'G', 'C']:

        # find which amino acid is encoded by a given codon
        def find_aa(codon_table,codon):
            for aa in codon_table.keys():
                if codon in codon_table[aa]:
                    return(aa)
            return(None)

        # count number of T's in codon
        def count_t(codon):
            count = 0
            for b in codon:
                if b == base:
                    count += 1
            return(count)

        # find most frequent codon, which encodes the same amino acid, but has less T's
        def find_replacement_codon(codon_table, aa):
            # sort by number of Ts (ascending) and codon frequency (descending)
            codons = sorted(codon_table[aa], key = lambda x: [count_t(x), -codon_table[aa][x]])
            # print(':'.join(['%s(%i,%f)' % (c,count_t(c),codon_table[aa][c]) for c in codons]))
            return(codons[0])

        # make sure that we get correct input
        seq_len = len(sequence)

        if seq_len % 3 != 0:
            print('[ERROR] Incorrect sequence length', seq_len)
            exit()

        # new sequence
        recoded_seq = ''

        for i in range(int(seq_len/3)):
            # current codon
            codon = sequence[i*3:i*3+3]

            # corresponding amino acid
            aa = find_aa(codon_table, codon)

            # make sure that codon table contains a given codon
            if aa == None:
                print('[ERROR] Could not find a matching codon in the codon_table', codon)
                exit()

            # find frequent codon with less T's (if available)
            new_codon = find_replacement_codon(codon_table, aa)

            # if found codon has the same number of T's, keep the original codon
            if count_t(codon) == count_t(new_codon):
                recoded_seq += codon
            else:
                # otherwise, use new codon
                recoded_seq += new_codon

            # # uncomment this line if you'd like to see the decision table
            #print(i, codon, aa, codon_table[aa][codon], new_codon, count_t(new_codon), codon_table[aa][new_codon], sep='\t')
                
        #make a new SeqRecord to process with Biopython
        new_bio_record = SeqRecord(Seq(recoded_seq),
            id=f"{seq_record.id}_{base}-depleted, GC%:{gc_fraction(recoded_seq):.2f}",  
            description="")
        record_list.append(new_bio_record)

    


#generate output fasta file for writing, and write
SeqIO.write(record_list, out_fasta_filename, "fasta")


#fuzzy search script is below
def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = bdi.ambiguous_dna_values
   r = []
   for i in product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r
ambiguous_motif = extend_ambiguous_dna(motif)


for count, record in enumerate(record_list,1):
    name, sequence = record.id, str(record.seq)    
    record_title =f"Record: {name}\n" 
    
    #screen printing
    print(count, record_title, end="")
    if count % 5 == 0:
        print("-----")

    new_record_title = f"{count}.{record_title}"   
    outmatch.write(new_record_title)    


    for subsequence in ambiguous_motif:
        title=(f"\tAmbiguous {subsequence}:\n")
        
        matches = find_near_matches(subsequence, sequence, 
            max_l_dist     = max_l_dist_, 
            max_deletions  = max_deletions_, 
            max_insertions = max_insertions_
        ) 
        

        #print(title, end="")
        outmatch.write(title)

        for item in matches:
            myline = (f"\t\t{item}\t{len(record)}\n") #prints to screen
            #print(myline, end="")
            outmatch.write(myline)
        if count % 5 == 0:
            outmatch.write("\n------------------\n\n")
        
outmatch.close()



    #------------------------------------------------------------


        # codon_table = {         # h_sapiens_9606
        #     '*': {'TAA': 0.3, 'TAG': 0.24, 'TGA': 0.47},
        #     'A': {'GCA': 0.23, 'GCC': 0.4, 'GCG': 0.11, 'GCT': 0.27},
        #     'C': {'TGC': 0.54, 'TGT': 0.46},
        #     'D': {'GAC': 0.54, 'GAT': 0.46},
        #     'E': {'GAA': 0.42, 'GAG': 0.58},
        #     'F': {'TTC': 0.54, 'TTT': 0.46},
        #     'G': {'GGA': 0.25, 'GGC': 0.34, 'GGG': 0.25, 'GGT': 0.16},
        #     'H': {'CAC': 0.58, 'CAT': 0.42},
        #     'I': {'ATA': 0.17, 'ATC': 0.47, 'ATT': 0.36},
        #     'K': {'AAA': 0.43, 'AAG': 0.57},
        #     'L': {'CTA': 0.07, 'CTC': 0.2, 'CTG': 0.4, 'CTT': 0.13, 'TTA': 0.08, 'TTG': 0.13},
        #     'M': {'ATG': 1.0},
        #     'N': {'AAC': 0.53, 'AAT': 0.47},
        #     'P': {'CCA': 0.28, 'CCC': 0.32, 'CCG': 0.11, 'CCT': 0.29},
        #     'Q': {'CAA': 0.27, 'CAG': 0.73},
        #     'R': {'AGA': 0.21, 'AGG': 0.21, 'CGA': 0.11, 'CGC': 0.18, 'CGG': 0.2, 'CGT': 0.08},
        #     'S': {'AGC': 0.24, 'AGT': 0.15, 'TCA': 0.15, 'TCC': 0.22, 'TCG': 0.05, 'TCT': 0.19},
        #     'T': {'ACA': 0.28, 'ACC': 0.36, 'ACG': 0.11, 'ACT': 0.25},
        #     'V': {'GTA': 0.12, 'GTC': 0.24, 'GTG': 0.46, 'GTT': 0.18},
        #     'W': {'TGG': 1.0},
        #     'Y': {'TAC': 0.56, 'TAT': 0.44}
        #     }
