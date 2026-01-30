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

