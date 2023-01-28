from Bio import SeqIO
from collections import Counter
'''
    Reads each tRNA sequence from fasta file
    Makes sliding window of based on window length entered by user.
'''

window_length = 4


def sliding_window(sequence, window_size):
    my_windows =[]
    if len(sequence) <= window_size:
        return my_windows.append(sequence)       
    
    for i in range(len(sequence)):
        my_windows.append(sequence[i:i+window_size])
    return my_windows


my_dict_list=[]
with open("hg38-tRNAs_2lines.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        windows = sliding_window(str(record.seq), window_length)
        kmers= Counter(windows)
        my_dict_list.append(kmers) #make a list of Counter dictionaries.


_count = Counter() #make empty counter
for item in my_dict_list:
    _count= _count+item

for k, v in _count.most_common():
    print(k, v)
