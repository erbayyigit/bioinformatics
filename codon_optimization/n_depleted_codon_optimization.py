"""
RNA Codon Optimizer - Nucleotide Depletion Tool

This script optimizes DNA/RNA sequences by replacing codons with synonymous alternatives
that have fewer instances of a target nucleotide, while maintaining amino acid sequences.

Features:
  - Automatic DNA to RNA conversion
  - Supports depletion of A, U, G, or C nucleotides
  - Calculates GC% before and after optimization
  - Uses E. coli codon usage frequencies for optimization
  - Works with both DNA and RNA FASTA files

Usage:
  python n_depleted_codon_optimization.py -i <input_file> -n <nucleotide>

Arguments:
  -i, --input        Input FASTA file (DNA or RNA sequence) [REQUIRED]
  -n, --nucleotide   Nucleotide to deplete (A, C, G, or U) [REQUIRED]
                     (Case-insensitive: a, u, g, c also accepted)

Examples:
  python n_depleted_codon_optimization.py -i sequences.fasta -n u
  python n_depleted_codon_optimization.py -i genes.fa -n G
  python n_depleted_codon_optimization.py -i transcript.fasta -n A

Output:
  FASTA format with header containing:
    - Sequence ID with '_<nucleotide>.depleted' suffix
    - GC% before optimization (GC_before=XX.XX%)
    - GC% after optimization (GC_after=XX.XX%)
"""

import argparse
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Optimize codons by depleting a specific nucleotide (RNA/U-based depletion)')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='Input FASTA file (DNA or RNA sequence)')
parser.add_argument('-n', '--nucleotide', type=str.upper, choices=['A', 'C', 'G', 'U'], 
                    required=True, help='Nucleotide to deplete in RNA (A, C, G, or U). Required.')
args = parser.parse_args()
input_file = args.input
target_nucleotide = args.nucleotide

# Helper functions to detect sequence type
def is_dna(seq):
    """Check if sequence is DNA (contains T) or RNA (contains U)"""
    return 'T' in seq and 'U' not in seq

def is_rna(seq):
    """Check if sequence is RNA (contains U) or DNA (contains T)"""
    return 'U' in seq and 'T' not in seq

# Read the first sequence to determine file type
records = list(SeqIO.parse(input_file, "fasta"))
if not records:
    print("[ERROR] FASTA file is empty.")
    exit()

first_seq = str(records[0].seq).upper()
if is_dna(first_seq):
    is_dna_file = True
elif is_rna(first_seq):
    is_dna_file = False
else:
    print("[ERROR] First sequence contains both T and U or neither. Cannot determine if DNA or RNA file.")
    exit()

# Process all sequences
for seq_record in records:
    sequence = str(seq_record.seq).upper()
    
    # Convert DNA to RNA if file is DNA
    if is_dna_file:
        sequence = sequence.replace('T', 'U')

    # e_coli_316407 (converted to RNA codons with U instead of T)
    ecoli_codon_table = {
        '*': {'UAA': 0.64, 'UAG': 0.07, 'UGA': 0.29},
        'A': {'GCA': 0.21, 'GCC': 0.27, 'GCG': 0.36, 'GCU': 0.16},
        'C': {'UGC': 0.56, 'UGU': 0.44},
        'D': {'GAC': 0.37, 'GAU': 0.63},
        'E': {'GAA': 0.69, 'GAG': 0.31},
        'F': {'UUC': 0.43, 'UUU': 0.57},
        'G': {'GGA': 0.11, 'GGC': 0.41, 'GGG': 0.15, 'GGU': 0.34},
        'H': {'CAC': 0.43, 'CAU': 0.57},
        'I': {'AUA': 0.07, 'AUC': 0.42, 'AUU': 0.51},
        'K': {'AAA': 0.76, 'AAG': 0.24},
        'L': {'CUA': 0.04, 'CUC': 0.1, 'CUG': 0.5, 'CUU': 0.1, 'UUA': 0.13, 'UUG': 0.13},
        'M': {'AUG': 1.0},
        'N': {'AAC': 0.55, 'AAU': 0.45},
        'P': {'CCA': 0.19, 'CCC': 0.12, 'CCG': 0.53, 'CCU': 0.16},
        'Q': {'CAA': 0.35, 'CAG': 0.65},
        'R': {'AGA': 0.04, 'AGG': 0.02, 'CGA': 0.06, 'CGC': 0.4, 'CGG': 0.1, 'CGU': 0.38},
        'S': {'AGC': 0.28, 'AGU': 0.15, 'UCA': 0.12, 'UCC': 0.15, 'UCG': 0.15, 'UCU': 0.15},
        'T': {'ACA': 0.13, 'ACC': 0.44, 'ACG': 0.27, 'ACU': 0.16},
        'V': {'GUA': 0.15, 'GUC': 0.22, 'GUG': 0.37, 'GUU': 0.26},
        'W': {'UGG': 1.0},
        'Y': {'UAC': 0.43, 'UAU': 0.57}
    }


    # find which amino acid is encoded by a given codon
    def find_aa(ecoli_codon_table,codon):
        for aa in ecoli_codon_table.keys():
            if codon in ecoli_codon_table[aa]:
                return(aa)
        return(None)

    # count occurrences of a specific nucleotide in codon
    def count_nucleotide(codon, nucleotide):
        count = 0
        for b in codon:
            if b == nucleotide:
                count += 1
        return(count)

    # find most frequent codon, which encodes the same amino acid, but has less of the target nucleotide
    def find_replacement_codon(ecoli_codon_table, aa, nucleotide):
        # sort by number of target nucleotides (ascending) and codon frequency (descending)
        codons = sorted(ecoli_codon_table[aa], key = lambda x: [count_nucleotide(x, nucleotide), -ecoli_codon_table[aa][x]])
        # print(':'.join(['%s(%i,%f)' % (c,count_nucleotide(c, nucleotide),ecoli_codon_table[aa][c]) for c in codons]))
        return(codons[0])

    # calculate GC% for a sequence
    def calculate_gc_percent(seq):
        gc_count = seq.count('G') + seq.count('C')
        return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

    # make sure that we get correct input
    seq_len = len(sequence)

    if seq_len % 3 != 0:
        print('[ERROR] Incorrect sequence length', seq_len)
        exit()

    # Calculate GC% before depletion
    gc_before = calculate_gc_percent(sequence)

    # new sequence
    recoded_seq = ''

    for i in range(int(seq_len/3)):
        # current codon
        codon = sequence[i*3:i*3+3]

        # corresponding amino acid
        aa = find_aa(ecoli_codon_table, codon)

        # make sure that codon table contains a given codon
        if aa == None:
            print('[ERROR] Could not find a matching codon in the codon_table', codon)
            exit()

        # find frequent codon with less target nucleotide (if available)
        new_codon = find_replacement_codon(ecoli_codon_table, aa, target_nucleotide)

        # if found codon has the same number of target nucleotides, keep the original codon
        if count_nucleotide(codon, target_nucleotide) == count_nucleotide(new_codon, target_nucleotide):
            recoded_seq += codon
        else:
            # otherwise, use new codon
            recoded_seq += new_codon

        # # uncomment this line if you'd like to see the decision table
        #print(i, codon, aa, ecoli_codon_table[aa][codon], new_codon, count_nucleotide(new_codon, target_nucleotide), ecoli_codon_table[aa][new_codon], sep='\t')
    
    # Calculate GC% after depletion
    gc_after = calculate_gc_percent(recoded_seq)
    
    print(f">{seq_record.id}_{target_nucleotide}.depleted| GC_before={gc_before:.2f}%| GC_after={gc_after:.2f}%")    
    print(recoded_seq)

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
