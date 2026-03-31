#!/usr/bin/env python3
"""
fasta-analyze-and-csv.py

Read protein records from a multirecord FASTA file,
compute basic protein properties, and write results as CSV.

Properties computed:
- Molecular Weight (kDa)
- Isoelectric point (pI)
- Charge at pH 7.0
- Cysteine count
- Amino acid sequence

Author: Erbay Yigit
Contact: erbayyigit@gmail.com
"""

from __future__ import annotations

import argparse
import csv
import sys
import os
import time
from typing import Sequence, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read protein records from a FASTA file and output analysis results as CSV."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to a FASTA file containing protein sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV file path. If omitted, a default based on input name will be used.",
    )
    return parser.parse_args(argv)


def read_fasta_records(fasta_path: str) -> List[SeqRecord]:
    """Read protein records from a FASTA file."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")
    
    return records


def compute_protein_metrics(record: SeqRecord) -> dict:
    """Compute protein metrics and return a dictionary for CSV output."""
    seq_str = str(record.seq)
    
    # Use record description as protein name, ID as fallback
    protein_name = (record.description or record.id).strip()
    
    try:
        P = ProteinAnalysis(seq_str)
        mw_kDa = P.molecular_weight() / 1000.0
        pI = P.isoelectric_point()
        cys = P.count_amino_acids().get("C", 0)
        charge_pH7 = P.charge_at_pH(7.0)
    except Exception as e:
        # Handle invalid sequences
        print(f"Warning: Could not analyze sequence {record.id}: {e}", file=sys.stderr)
        return {
            "sequence_id": record.id,
            "protein_name": protein_name,
            "length": len(seq_str),
            "mw_kDa": "N/A",
            "pI": "N/A",
            "cys": "N/A",
            "charge_pH7": "N/A",
            "protein": seq_str,
        }

    return {
        "sequence_id": record.id,
        "protein_name": protein_name,
        "length": len(seq_str),
        "mw_kDa": mw_kDa,
        "pI": pI,
        "cys": cys,
        "charge_pH7": charge_pH7,
        "protein": seq_str,
    }


def write_csv(records: Sequence[SeqRecord], out_path: str) -> None:
    """Write computed metrics to CSV."""
    fieldnames = [
        "sequence_id",
        "protein_name",
        "length",
        "mw_kDa",
        "pI",
        "charge_pH7",
        "cys",
        "protein",
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for record in records:
            m = compute_protein_metrics(record)
            # Format numeric values if they exist
            if isinstance(m["mw_kDa"], float):
                m["mw_kDa"] = f"{m['mw_kDa']:.2f}"
            if isinstance(m["pI"], float):
                m["pI"] = f"{m['pI']:.2f}"
            if isinstance(m["charge_pH7"], float):
                m["charge_pH7"] = f"{m['charge_pH7']:.2f}"
            writer.writerow(m)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)

    # Read FASTA file
    try:
        records = read_fasta_records(args.input)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine output path
    if args.output is None:
        input_name = os.path.splitext(os.path.basename(args.input))[0]
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_dir = os.path.join(os.path.dirname(__file__), "output")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        args.output = os.path.join(output_dir, f"{input_name}_analysis_{timestamp}.csv")
    else:
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

    # Write results
    write_csv(records, args.output)
    print(f"Analyzed {len(records)} sequences and wrote results to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
