import subprocess
from datetime import datetime
import sys

# Get input file from command line (like the $1 in bash)
input_file = sys.argv[1]


# Create the timestamped filename
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
base_name = input_file.rsplit('.', 1)[0]
output_file = f"results_{base_name}_{timestamp}.fasta"


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

















# #!/bin/bash

# # $1 represents the first argument you type after the script name
# INPUT_FILE=$1

# # This creates a timestamp like 2026-01-29_16-20
# TIMESTAMP=$(date +"%Y-%m-%d_%H-%M")

# # This removes the extension from input and builds: results_filename_timestamp.fasta
# BASE_NAME="${INPUT_FILE%.*}"
# OUT_FILE="results_${BASE_NAME}_${TIMESTAMP}.fasta"

# while read -r line <&3; do
#     # It is good practice to echo the ID so you can see progress in the terminal
#     echo "Fetching $line..."
#     efetch -db protein -id "$line" -format fasta
# done 3< "$INPUT_FILE" > "$OUT_FILE"

# echo "Process complete. Output saved to: $OUT_FILE"

