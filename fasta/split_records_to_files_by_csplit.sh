#yigit@neb.com 
#Dec 31, 2022

# This tutorial for GNU csplit not mac csplit
# Use it from Core Utilities as gcplit

# takes each fasta record and split into one recor per file
# nomatter how many records you have
#the first file will be emtpy
gcsplit -f chr temp.fasta "/^>/" {*}

#this will do make 6 files with single record
# and remaining 45 records out of 50 will be in the 7th file.
# remember the first line number start with 0.
gcsplit input.fasta "/^>/" {5}







