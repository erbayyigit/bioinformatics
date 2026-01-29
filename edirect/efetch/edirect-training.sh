
#Example-1
#Finde most prolific authors who published on m6A
esearch -db pubmed -query "m6A"|\
efetch -format xml|\
xtract -pattern PubmedArticle -block Author \
-sep " " -tab "\n" -element LastName,Initials|\
sort-uniq-count-rank
