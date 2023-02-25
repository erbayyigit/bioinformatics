'''
   It reads genbank IDs from a file
   fetches records from internet
   prints annotation etc to CSV file
   or
   generates a FASTA file

   yigit@neb.com
'''

from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

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


with open('resuls.csv', 'w', newline='') as csvfile:
   csv_writer = csv.writer(csvfile, delimiter = ',')
   csv_writer.writerow(['accession', 'protein', 'organism',  'mw', 'isoep', 'cys', 'charge@pH7'])
   
   for record in seq_records:
      organism =record.annotations['organism']      
      P = ProteinAnalysis(record.seq)
      mw = f'{(P.molecular_weight()/1000):.2f}'
      isoep = f'{P.isoelectric_point():.2f}'
      cys = P.count_amino_acids()['C']
      charge=f'{P.charge_at_pH(7.0):.2f}'   
      csv_writer.writerow([record.id, record.seq, organism,  mw, isoep, cys, charge])
      #use below to write fasta
      # print(f'>{record.id} {mw} {isoep} pH7:{charge} Cys:{cystesins}\n{record.seq}')
