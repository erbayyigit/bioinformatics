#!/usr/bin/env python3
import operator
import os
import sys
from Bio import Seq, SeqIO, Entrez

def fetch_genbank_batch(id_file, output_file):
    """Fetch FASTA records from GenBank (NCBI) for a list of IDs."""
    
    # Check if input file exists
    if not os.path.exists(id_file):
        print(f"Error: Input file '{id_file}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Check if output file already exists
    if os.path.exists(output_file):
        print(f"Error: Output file '{output_file}' already exists.", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # Set Entrez email for NCBI
    Entrez.email = "erbayyigit@gmail.com"
    
    # Read IDs from file
    ids = []
    successful = 0
    failed = 0
    
    print(f"Reading IDs from: {id_file}")
    with open(id_file) as fh:
        for line in fh:
            line = line.strip()
            if line:
                ids.extend(line.split())
    
    if not ids:
        print("Error: No IDs found in input file.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(ids)} IDs to retrieve\n")
    
    try:
        # Fetch records from GenBank
        print("Fetching sequences from GenBank...")
        handle = Entrez.efetch(db="protein", rettype="fasta", id=ids)
        seq_records = list(SeqIO.parse(handle, "fasta"))
        
        if not seq_records:
            print("Error: No sequences retrieved from GenBank.", file=sys.stderr)
            sys.exit(1)
        
        # Write to output file
        with open(output_file, 'w') as f_out:
            for seq_record in seq_records:
                SeqIO.write(seq_record, f_out, "fasta")
                successful += 1
        
        print(f"✓ Successfully retrieved {successful} sequences")
        
    except Exception as e:
        print(f"✗ Error fetching from GenBank: {e}", file=sys.stderr)
        failed = len(ids)
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    # Get input filename from command line argument
    if len(sys.argv) < 2:
        print("Usage: python retrieve_records_by_genebank_id.py <input_file> [output_file]")
        print("Example: python retrieve_records_by_genebank_id.py input/genbank_ids.txt output/results.fasta")
        sys.exit(1)
    
    id_file = sys.argv[1]
    
    # Get output filename or use default
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = "output/results.fasta"
    
    fetch_genbank_batch(id_file, output_file)

