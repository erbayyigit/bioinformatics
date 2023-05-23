#convert pandas dataframe to FASTA file
import pandas as pd

df = pd.read_csv('yL_syngen_database.csv', usecols = ['syn_gene_id','seq', 
    'promoter'])

#select columns that has t7 promoter
df=df.query("promoter == 'T7 promoter'")
#drop Ids from sequences
df=df.dropna()
df=df.drop(['promoter'], axis=1)

#converts things to list
record_id =  df['syn_gene_id'].tolist()
record_id = [">"+element for element in record_id]
record_seq = df['seq'].tolist()

my_dictionary = dict()
my_dictionary=dict(zip((record_id), record_seq))

for key in my_dictionary:
    print(key+"\n"+ my_dictionary[key])

