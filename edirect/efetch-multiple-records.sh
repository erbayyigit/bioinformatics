while read line <&3; do
	efetch -db protein -id "$line" -format fasta
done 3< sorted_thermophilic_ids.txt > thermophilic_mazf_all.fasta
