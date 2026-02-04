import operator
from Bio import Seq, SeqIO, Entrez
#this is just to get gb ids into a list
with open('rosalind_frmt.txt') as fh:
    for line in fh:
        line = line.strip()
        ids=line.split()

#this fetches a list of genebank ids in fasta format.
Entrez.email="erbayyigit@gmail.com"
handle=Entrez.efetch(
    db="nucleotide", rettype="fasta", 
    id=ids) 
seq_records = list(SeqIO.parse(handle, "fasta"))

#I do thi so that I can later sort dic keys by lenght
my_dict={}
for seq_record in seq_records:
    my_dict[len(seq_record.seq)] = f'>{seq_record.description}\n{seq_record.seq}'

#my goal is to choose shortest fasta record from given ones.
sorted_dic = sorted(my_dict.items(), key=operator.itemgetter(0))
print(sorted_dic[0][1])

