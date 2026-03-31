#!/usr/bin/env python3
import requests
from Bio import SeqIO
from io import StringIO
import time
import os
import sys

def fetch_uniprot_batch(id_file, output_file):
    """Fetch FASTA records from UniProt for a list of IDs."""
    
    # Check if input file exists
    if not os.path.exists(id_file):
        print(f"Error: Input file '{id_file}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    successful = 0
    failed = 0
    
    with open(id_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line_num, line in enumerate(f_in, 1):
            accession = line.strip()
            if not accession:
                continue
            
            # Construct the direct UniProt URL
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
            
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    # Use StringIO to let Biopython "read" the text string as a file
                    fasta_data = StringIO(response.text)
                    record = SeqIO.read(fasta_data, "fasta")
                    
                    # Write to your output file
                    SeqIO.write(record, f_out, "fasta")
                    print(f"✓ Line {line_num}: Successfully retrieved {accession}")
                    successful += 1
                else:
                    print(f"✗ Line {line_num}: Failed to find {accession} (Status: {response.status_code})", file=sys.stderr)
                    failed += 1
                
                # Polite "sleep" to avoid hitting the API too fast
                time.sleep(0.1) 
                
            except Exception as e:
                print(f"✗ Line {line_num}: Error processing {accession}: {e}", file=sys.stderr)
                failed += 1
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    fetch_uniprot_batch("input/uniprot_ids.txt", "output/results.fasta")
