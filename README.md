# edirect-utils

A collection of Python utilities for working with NCBI EDirect and protein sequence data.

## Contents

- **efetch-proteins-by-id.py**  
  Fetches protein sequences from NCBI using a list of accession IDs and saves them as a FASTA file.

- **efetch-proteins-by-id-and-analyze.py**  
  Fetches protein records from NCBI, computes basic protein properties (e.g., molecular weight, pI, cysteine count), and outputs results as CSV or annotated FASTA.

- **(future scripts...)**  
  This repository will expand with more bioinformatics utilities.

## Requirements

- Python 3.x (created on Python 3.10)
- [NCBI EDirect tools](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (for efetch-proteins-by-id.py)
- [Biopython](https://biopython.org/) (for analysis scripts)

## Usage

Each script contains usage instructions in its header and can be run from the command line. Example:

```bash
python efetch-proteins-by-id.py input_ids.txt
python efetch-proteins-by-id-and-analyze.py -i input_ids.txt -f CSV
```

## Output

Results are saved in the `/output` directory, with filenames including timestamps to prevent overwriting.

