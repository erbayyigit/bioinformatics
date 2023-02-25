from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

temp_list=[] # get list of genbank IDs srater with WP_
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
   #print(record.id,  record.annotations['organism'], record.seq)
   P = ProteinAnalysis(record.seq)
   mw = f'{(P.molecular_weight()/1000):.2f}kDa'
   isoep = f'{P.isoelectric_point():.2f}'
   cystesins = P.count_amino_acids()['C']
   charge=f'{P.charge_at_pH(7.0):.2f}'
   
   print(f'>{record.id} {mw} {isoep} pH7:{charge} Cys:{cystesins}\n{record.seq}')

