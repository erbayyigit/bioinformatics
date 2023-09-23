while read line <&3; do
	curl -s "https://rest.uniprot.org/uniprotkb/$acc.fasta"
done 3 < uniprot_ids.txt > outputfile.fasta

