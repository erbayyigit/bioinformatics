import os
#write each record to a new file from a multirecord fasta file

with open("input.fasta") as f:
	for line in f:
		if line.startswith('>'):
			outfilename = f"{(line.rstrip())}.fasta"
			if not os.path.exists(outfilename):
				out = open(outfilename, 'w+')		
			else:
				print('File already exists!')
		out.write(line)
