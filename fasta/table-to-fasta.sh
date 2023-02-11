#by Erbay Yigit
#convert two column tabular file to FASTA
awk -F"," '{print ">"$1"\n"$2}' file.csv > outputFile.fasta
#input
#seqID, AGTCG...

#output format
#>seqID
#AGTCG...


