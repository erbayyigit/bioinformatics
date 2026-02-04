"""
Author: Erbay Yigit
Email: erbayyigit@gmail.com

Purpose:
--------
This script retrieves protein sequences from NCBI using a list of protein
accession IDs and saves them as a single FASTA file.

It reads an input text file containing one protein ID per line, fetches
each corresponding sequence via NCBI's `efetch` command-line utility,
and writes all retrieved sequences to a timestamped FASTA output file.

Usage:
------
python fetch_proteins.py <input_file>

Input:
------
- <input_file>: A text file containing one protein accession ID per line.

Output:
-------
- A FASTA file named:
  result_<input_filename>_<YYYYMMDD_HHMM>.fasta

Requirements:
-------------
- NCBI EDirect tools (`efetch`) must be installed and accessible in PATH.
- Python 3.x

Notes:
------
- Empty lines in the input file are ignored.
- The timestamp prevents accidental overwriting of previous results.
"""


import os, subprocess, sys
from datetime import datetime


# Get input file from command line
input_file = sys.argv[1]

# 1. Get just the filename (e.g., "input.txt" from "samples/input.txt")
file_only = os.path.basename(input_file)

# 2. Get the name without the extension (e.g., "input" from "input.txt")
base_name = os.path.splitext(file_only)[0]
# Create the timestamped filename

timestamp = datetime.now().strftime("%Y%m%d_%H%M")

output_file = f"result_{base_name}_{timestamp}.fasta"


with open(input_file, 'r') as ids, open(output_file, 'w') as out:
    for line in ids:
        protein_id = line.strip()
        if protein_id:
            print(f"Fetching {protein_id}...")
            # This runs the actual efetch command
            result = subprocess.run(
                ["efetch", "-db", "protein", "-id", protein_id, "-format", "fasta"],
                capture_output=True,
                text=True
            )
            out.write(result.stdout)

print(f"Success! Saved to {output_file}")

