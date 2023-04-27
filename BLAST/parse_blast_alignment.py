'''
    This script parses BLAST alignment XML files.
    It greps the "DEscription" field.
    Prints unique hits
    yigit@neb.com 2023.04.19
'''
from Bio.Blast import NCBIXML
from Bio import SearchIO


result_handle = open("alignment.xml")
blast_records = list(NCBIXML.parse(result_handle))

my_description_list =[]
for blast_record in blast_records:
    for one_hit in blast_record.alignments:
        #print(one_hit.hit_id)
        myhit = one_hit.hit_def
        description = myhit.split(' [', 1) #split before species name at first occurence
        my_description_list.append(description[0])
    my_set_of_descriptons = set(my_description_list)
    joined_descriptions = '; '.join(my_set_of_descriptons)
    joined_for_csv = blast_record.query+"\t"+joined_descriptions
    # print(blast_record.query, ",", joined_descriptions)
    print(joined_for_csv)
    my_description_list=[] #reset the

