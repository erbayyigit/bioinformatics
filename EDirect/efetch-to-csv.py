'''
   It reads genbank IDs from a file
   fetches records from internet
   prints annotation etc to CSV file
   or
   generates a FASTA file
   yigit@neb.com
   BioPython Version 1.81
'''
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
'''
# use the following code if
# gbk ideas are enumated.
temp_list=[]
with open('gene-list.txt', 'r') as infile:
   for gb_id in infile:
      gb_id =gb_id.strip()
      #print(gb_id)
      _, gb_id = gb_id.split("_", 1)
      temp_list.append(gb_id)
'''
MBP = ProteinAnalysis(
"KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDA"
"VRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLT"
"FLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGL"
"EAVNKDKPLGAVALKSYEEELVKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTGSGSGSLVPRGSSGSSHHHHHH")
mbp_mw    = MBP.molecular_weight()
mbp_isoep = MBP.isoelectric_point()

gbk_id_list = []
with open('typ3a-gbk-ids.txt', 'r') as infile:
   for gb_id in infile:
      gb_id =gb_id.strip()
      gbk_id_list.append(gb_id)

list_of_WP_ids = set(gbk_id_list) 
Entrez.email = "erbayyigit@gmail.com"
handle = Entrez.efetch(db="protein", rettype="gb",  retmode="text", id=list_of_WP_ids)
seq_records = list(SeqIO.parse(handle, "gb"))

with open('results.csv', 'w', newline='') as csvfile:
   csv_writer = csv.writer(csvfile, delimiter = ',')
   csv_writer.writerow(['accession', 'protein', 'organism', 'mw', 'mbp_fusion','isoep', 'cys', 'charge@pH7'])
   
   for record in seq_records:
      organism =record.annotations['organism']      
      P = ProteinAnalysis(record.seq)
      mw = f'{(P.molecular_weight()/1000):.2f}'
      mbp_fusion = float(mw) + float(42.42)
      isoep = f'{P.isoelectric_point():.2f}'
      cys = P.count_amino_acids()['C']
      charge=f'{P.charge_at_pH(7.0):.2f}'   
      csv_writer.writerow([record.id, record.seq, organism, mw, mbp_fusion, isoep, cys, charge])
      #use below to write fasta
      # print(f'>{record.id} {mw} {isoep} pH7:{charge} Cys:{cystesins}\n{record.seq}')

