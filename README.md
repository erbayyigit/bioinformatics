# Bioinformatics Sequence Retrieval Scripts

A collection of Python utilities for fetching protein and nucleotide sequences from public databases (UniProt and GenBank/NCBI).

## Scripts

### 1. retrieve_records_by_uniprot_ids.py

Fetches protein sequences from UniProt REST API using a list of UniProt accession IDs.

**Features:**
- Reads UniProt IDs from input file
- Fetches FASTA records from the UniProt REST API
- Writes results to output FASTA file
- Automatic output directory creation
- Progress tracking and summary statistics
- Error handling with detailed logging

**Usage:**
```bash
python scripts/retrieve_records_by_uniprot_ids.py
```

**Input:** `input/uniprot_ids.txt` (one ID per line)  
**Output:** `output/results.fasta`

**Example input file:**
```
P12345
A0A024B7W1
Q9Y5K6
```

---

### 2. retrieve_records_by_genebank_id.py

Fetches protein sequences from GenBank (NCBI) using a list of GenBank accession IDs.

**Features:**
- Reads GenBank IDs from command-line argument
- Fetches FASTA records from NCBI Entrez API
- Writes results to output FASTA file (default: `output/results.fasta`)
- Automatic output directory creation
- Prevents accidental file overwriting
- Progress tracking and summary statistics
- Error handling with detailed logging

**Usage:**
```bash
python scripts/retrieve_records_by_genebank_id.py <input_file> [output_file]
```

**Examples:**
```bash
# Use default output location
python scripts/retrieve_records_by_genebank_id.py input/genbank_ids.txt

# Specify custom output location
python scripts/retrieve_records_by_genebank_id.py input/genbank_ids.txt output/custom_results.fasta
```

**Input:** Comma or space-separated GenBank IDs (one per line or multiple per line)  
**Output:** `output/results.fasta` (or custom path)

**Example input file:**
```
NM_000001
NC_000001
AF123456
```

---

## Project Structure

```
bioinformatics/
├── scripts/
│   ├── retrieve_records_by_uniprot_ids.py
│   └── retrieve_records_by_genebank_id.py
├── input/
│   ├── uniprot_ids.txt
│   └── genbank_ids.txt
├── output/
│   └── results.fasta
└── README.md
```

---

## Requirements

- Python 3.6+
- Dependencies:
  - `requests` - For HTTP requests
  - `biopython` - For sequence parsing and writing
  - `Bio.Entrez` - For NCBI database access (included in biopython)

**Install dependencies:**
```bash
pip install requests biopython
```

---

## Configuration

### UniProt Script
The UniProt script uses hardcoded paths:
- Input: `input/uniprot_ids.txt`
- Output: `output/results.fasta`

### GenBank Script
The GenBank script requires command-line arguments:
```bash
python scripts/retrieve_records_by_genebank_id.py <input_file> [output_file]
```

The NCBI Entrez email is set to: `erbayyigit@gmail.com`  
(Modify in the script if needed)

---

## Features

✅ Automatic directory creation  
✅ Error handling and validation  
✅ Input file existence checking  
✅ Duplicate file prevention (GenBank script)  
✅ Progress tracking during retrieval  
✅ Summary statistics after completion  
✅ Detailed logging to stderr for errors  

---

## Error Handling

Both scripts include:
- **Input validation:** Checks if input file exists
- **Output validation:** Prevents overwriting existing files (GenBank)
- **API error handling:** Catches HTTP errors and parsing exceptions
- **Polite rate limiting:** 0.1 second delay between API requests

---

## Output Format

All output is in FASTA format:
```
>sequence_id description
MKVLWAALLV...
```

Each sequence retrieved is appended to the output file.

---

## Notes

- The UniProt script uses fixed file paths. Modify the script if you need different paths.
- The GenBank script accepts input filename as a command-line argument.
- Both scripts create the `output/` directory automatically if it doesn't exist.
- Rate limiting is set to 0.1 seconds between requests to be respectful to public APIs.
- An NCBI Entrez email is required for GenBank queries. Update it in the script if needed.

---

## License

This project is open source.

