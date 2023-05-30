import itertools
'''
	awk 1 ORS='' complete_5_25_1_100.txt |wc
	improvement needed:
	make sure this works with FASTA file
	yigit@neb.com May 29, 2023
	Boxford
'''

size_of_kmer = 7


#***********************************************************************************
# let's get curlcake output file and concatanate to single string thus you
# can search for pentamers in it if necessary. This will slurp the content into a single string
with open('../data/pey43_5mers_IVT.txt') as x:
	curlcakeRNA = x.read().replace('\n','')	# while slurping remove new line characters as well.
#***********************************************************************************


#***********************************************************************************
#Make kmers from sequence provided from file
def split_to_kmers(dna, kmer_size):
	kmers=[]
	for start in range(0, len(dna) - kmer_size + 1, 1):
		kmer = dna[start:start+kmer_size]
		kmers.append(kmer)
	return kmers

kmers_of_mysequence = split_to_kmers(curlcakeRNA, size_of_kmer)
uniq_kmers_of_mysequence =set(kmers_of_mysequence)
#***********************************************************************************



#***********************************************************************************
# HERE WE ARE MAKING ALL THEORITICAL PERMUTATION OF PENTAMERS
x = 'ACTG' # Do not change this. Itertools use it to make permutation, RNA bases.
theoretical_kmers = []
for output in itertools.product(x, repeat = size_of_kmer): # create pentamers using A,C,U,G
    theoretical_kmers.append(''.join(output))
#print(theoretical_kmers)
#print(len(theoretical_kmers))
#***********************************************************************************

# YOU CAN CHECK IF YOU HAVE ALL THE PENTAMERS IN CURLCAKE SEQUENCE.
for kmer in uniq_kmers_of_mysequence:
	if kmer in theoretical_kmers:
		print(kmer)



