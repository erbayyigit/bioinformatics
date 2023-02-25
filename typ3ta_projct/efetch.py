'''
It will fetch protein record based on genbank IDs.
List of attributes can be found here.
https://biopython.org/docs/1.74/api/Bio.GenBank.Record.html
yigit@neb.com 
'''

from Bio import SeqIO, Entrez

temp_list=[]
with open('type3tas_accession_numbers.txt', 'r') as infile:
   for gb_id in infile:
      gb_id =gb_id.strip()
      _, gb_id = gb_id.split("_",1)
      temp_list.append(gb_id)

list_of_WP_ids = set(temp_list)      
Entrez.email = "erbayyigit@gmail.com"

handle = Entrez.efetch(db="protein", rettype="gb",  retmode="text", id=list_of_WP_ids)
seq_records = list(SeqIO.parse(handle, "gb"))
for record in seq_records:
   print(record.id,  record.annotations['organism'], record.seq)

