#yigit_python_library
# January 27, 2023
#USAGE
# 	import sys
#	sys.path.insert(0, '/Users/supertech/Dropbox/comp-python')
# 	from yigit_python_lib import number_u
# 	print(number_u("AUUGC"))
#----------------------------------------------------------------------------
def base_frequency(infile:str, format:str):
	''' 
	Script: This script will pring base frequencey in the fasta record title.
	Read gbk or fasta file.
	Argument Usage:
	infile= myfile.gbk
	format = genebank or format=fasta
	Call example: base_frequency("AB610939-AB610950.gb","genbank")
	'''
	from Bio import SeqIO
	import re
	import textwrap
	for index, record in enumerate(SeqIO.parse(infile, format)): 
	    dna = str(record.seq)
	    g_count = dna.count('G')
	    c_count = dna.count('C')
	    a_count = dna.count('A')
	    rna = dna.replace('T', 'U' )
	    u_count = rna.count('U')
	    
	    g_pct= g_count/len(dna)
	    c_pct= c_count/len(dna)
	    a_pct= a_count/len(dna)
	    u_pct= u_count/len(dna)

	    seqlength = len(dna)
	    fasta_record=textwrap.fill(dna, width=60)

	    res=(f">{record.id}_{seqlength}nt g:{g_pct:.2f} c:{c_pct:.2f} \
	    	a:{a_pct:.2f} u:{u_pct:.2f}\n{fasta_record}")
	    return res


def double_stranded(myseq):
	'''
	myseq should come from fasta file. you can use read_fasta function for this.
	'''
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	
	sense = "Sense:\t" + myseq.seq
	count = len(myseq.seq)
	align_pipe = "\t" + '|' * count
	comp = "Comp:\t" + myseq.seq.complement()
	return sense +"\n"+ align_pipe +"\n" + comp


# def revcom(myseq):
# 	from Bio.Seq import Seq
# 	from Bio.Alphabet import IUPAC

def number_u(dna):
	'''
	countnumber "U"
	'''
	import re	# required for regex
	dna = dna.upper()
	rna = dna.replace('T', 'U' )
	u_count = rna.count('U')
	rnalength = len(rna)
	pct_u = 100*(u_count/rnalength)
	result = "Number of Us\t:" + str(u_count)+ "\n"+					\
	"RNA legnth\t:"+str(rnalength) + "\n" + \
	"Pct U\t:"+ str(pct_u) 
	return result


#this function calculate gc percent of DNA
def pct_gc(dna, sig_figs=2):
	length = len(dna)
	g_count = dna.upper().count('G')
	c_count = dna.upper().count('C')
	gc_content = (g_count+c_count)/length
	return "pct_gc\t:" + str(round(gc_content, sig_figs))



#------********** DRAFT
# def read_fasta_file(infile, format="fasta"):
# 	''' 
# 		This is going to read fasta file many records

# 	'''
# 	from Bio import SeqIO

#------********** DRAFT


def read_single_fasta_record(infile, format="fasta"):
	''' 
	read file containing **single fasta entry**. 
	'''
	from Bio import SeqIO
	record = SeqIO.read(infile, format)
	return record


def read_single_fasta_rec_as_string(infile, format="fasta"):
	''' 
	read file containing **single fasta entry**. 
	usage example:
	record = read_single_fasta(infile="xp12.fasta")
	print(record)
	return plain string !!
	'''
	from Bio import SeqIO
	record = SeqIO.read(infile, format)
	plain_string = str(record.seq)  #this converts sequence object to plain string
	return plain_string  #now you have plain dna or rna string do whatever you want with it.
#----------------------------------------------------------------------------


def reverse_complement_rna(rna):
    ''' 
		get reverse coplement of an RNA sequence
    '''
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return ''.join([complement[base] for base in rna[::-1]])


def sliding_window(sequence, window_size):
    '''
    You can use this function to make sliding windows from fasta
    sequence: DNA or RNA sequence. It can come from biopyton record.seq
    window_size: it's the size of window, such as 10nt
    my_windows contains list of kmers
	USAGE:
	windows = sliding_window(str(record.seq), window_length)
	kmers= Counter(windows)
	my_dict_list=[]
	my_dict_list.append(kmers)
    '''
    my_windows =[]
    if len(sequence) <= window_size:
        return my_windows.append(sequence)       
    
    for i in range(len(sequence)):
        my_windows.append(sequence[i:i+window_size])
    return my_windows


	
