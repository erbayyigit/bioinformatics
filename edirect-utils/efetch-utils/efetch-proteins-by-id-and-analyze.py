#!/usr/bin/env python3
"""
efetch-and-analyze.py

Fetch protein records from NCBI (Entrez) using a list of accession/GenBank IDs,
compute basic protein properties, and write results as either CSV or FASTA.

Adds:
- protein_name (from record.description and/or GB qualifiers)
- gene_name (from GenBank feature qualifiers when present)

Author: Erbay Yigit
Contact: erbayyigit@gmail.com
"""

from __future__ import annotations

import argparse
import csv
import sys
import os
import time
from typing import Iterable, List, Sequence, Tuple

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch NCBI protein records from a list of IDs and output CSV or FASTA."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to a text file containing one protein ID per line.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["CSV", "FASTA"],
        default="CSV",
        help="Output format (default: CSV).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output file path. If omitted, a default will be chosen based on format.",
    )
    parser.add_argument(
        "--email",
        default="erbayyigit@gmail.com",
        help="Email address required by NCBI Entrez.",
    )
    parser.add_argument(
        "--db",
        default="protein",
        help="NCBI Entrez database to query (default: protein).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=200,
        help="How many IDs to fetch per Entrez request (default: 200).",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.34,
        help="Seconds to sleep between chunks (default: 0.34).",
    )
    parser.add_argument(
        "--fasta-annotate",
        action="store_true",
        help="FASTA headers include computed properties (mw_kDa, pI, charge@pH7, Cys) plus names.",
    )
    parser.add_argument(
        "--fasta-name",
        action="store_true",
        help="Include gene/protein names in FASTA header (even if --fasta-annotate is not set).",
    )
    return parser.parse_args(argv)


def read_id_file(path: str) -> List[str]:
    """Read IDs from a file, one per line; ignore blanks; de-duplicate preserving order."""
    seen = set()
    ids: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            pid = line.strip()
            if not pid:
                continue
            if pid not in seen:
                seen.add(pid)
                ids.append(pid)
    if not ids:
        raise ValueError(f"No IDs found in input file: {path}")
    return ids


def chunked(items: Sequence[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield list(items[i : i + size])


def fetch_gb_records(
    protein_ids: Sequence[str],
    db: str,
    email: str,
    chunk_size: int = 200,
    sleep_s: float = 0.34,
) -> List[SeqRecord]:
    """Fetch GenBank records from NCBI and parse into SeqRecord objects."""
    Entrez.email = email
    all_records: List[SeqRecord] = []

    for id_chunk in chunked(protein_ids, chunk_size):
        handle = Entrez.efetch(db=db, rettype="gb", retmode="text", id=",".join(id_chunk))
        try:
            all_records.extend(list(SeqIO.parse(handle, "gb")))
        finally:
            handle.close()

        if sleep_s and len(protein_ids) > chunk_size:
            time.sleep(sleep_s)

    return all_records


def _first_qualifier(feature, keys: Sequence[str]) -> str:
    """Return the first found qualifier value for any key in keys, else ''."""
    quals = getattr(feature, "qualifiers", {}) or {}
    for k in keys:
        if k in quals and quals[k]:
            v = quals[k][0]
            return v.strip()
    return ""


def extract_names(record: SeqRecord) -> Tuple[str, str]:
    """
    Extract (protein_name, gene_name) from a GenBank protein SeqRecord.

    Notes:
    - For protein records, gene info is often in the CDS feature qualifiers (gene/locus_tag).
    - Protein name is often in record.description and/or CDS 'product' qualifier.
    """
    protein_name = (record.description or "").strip()
    gene_name = ""

    # Try features (common for GenBank/INSDC records)
    for feat in getattr(record, "features", []) or []:
        if feat.type == "CDS":
            # Protein name from product if it looks better than description
            product = _first_qualifier(feat, ["product"])
            if product and (not protein_name or protein_name == record.id):
                protein_name = product

            # Gene name preference: gene, then locus_tag
            gene_name = _first_qualifier(feat, ["gene", "locus_tag"])
            break

    # If description is overly generic or empty, keep whatever we found
    if protein_name == record.id:
        # Some records have id repeated as description
        protein_name = ""

    return protein_name, gene_name


def compute_protein_metrics(record: SeqRecord) -> dict:
    """Compute protein metrics and return a dictionary for CSV/FASTA annotation."""
    organism = record.annotations.get("organism", "") or ""
    seq_str = str(record.seq)

    protein_name, gene_name = extract_names(record)

    P = ProteinAnalysis(seq_str)
    mw_kDa = P.molecular_weight() / 1000.0
    pI = P.isoelectric_point()
    cys = P.count_amino_acids().get("C", 0)
    charge_pH7 = P.charge_at_pH(7.0)

    return {
        "accession": record.id,
        "protein": seq_str,
        "organism": organism,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "mw_kDa": mw_kDa,
        "pI": pI,
        "cys": cys,
        "charge_pH7": charge_pH7,
    }


def write_csv(records: Sequence[SeqRecord], out_path: str) -> None:
    """Write computed metrics to CSV, including gene/protein names."""
    fieldnames = [
        "accession",
        "protein_name",
        "gene_name",
        "organism",
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
            m["mw_kDa"] = f"{m['mw_kDa']:.2f}"
            m["pI"] = f"{m['pI']:.2f}"
            m["charge_pH7"] = f"{m['charge_pH7']:.2f}"
            writer.writerow(m)


def write_fasta(
    records: Sequence[SeqRecord],
    out_path: str,
    annotate: bool = False,
    include_names: bool = False,
) -> None:
    """
    Write sequences to FASTA.

    - annotate=True: add computed properties to header (and names)
    - include_names=True: add gene/protein names to header even without full annotation
    """
    with open(out_path, "w", encoding="utf-8") as out:
        for record in records:
            if annotate or include_names:
                m = compute_protein_metrics(record)
                protein_name = m["protein_name"].replace(" ", "_") if m["protein_name"] else ""
                gene_name = m["gene_name"].replace(" ", "_") if m["gene_name"] else ""

                name_bits = []
                if protein_name:
                    name_bits.append(f"protein:{protein_name}")
                if gene_name:
                    name_bits.append(f"gene:{gene_name}")

                if annotate:
                    header = (
                        f"{m['accession']} "
                        + " ".join(name_bits + [
                            f"mw_kDa:{m['mw_kDa']:.2f}",
                            f"pI:{m['pI']:.2f}",
                            f"pH7:{m['charge_pH7']:.2f}",
                            f"Cys:{m['cys']}",
                            f"organism:{m['organism'].replace(' ', '_')}",
                        ])
                    ).strip()
                else:
                    # Only names (no computed properties)
                    header = f"{record.id} " + " ".join(name_bits)
                    header = header.strip()
            else:
                header = record.id

            out.write(f">{header}\n{record.seq}\n")


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)


    # Always write to the /output directory relative to this script
    output_dir = os.path.join(os.path.dirname(__file__), "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    if args.output is None:
        base = "results.csv" if args.format == "CSV" else "results.fasta"
    else:
        base = os.path.basename(args.output)
    name, ext = os.path.splitext(base)
    out_filename = f"{name}_{timestamp}{ext}"
    args.output = os.path.join(output_dir, out_filename)

    protein_ids = read_id_file(args.input)
    records = fetch_gb_records(
        protein_ids=protein_ids,
        db=args.db,
        email=args.email,
        chunk_size=args.chunk_size,
        sleep_s=args.sleep,
    )

    if not records:
        print("No records returned. Check IDs and database.", file=sys.stderr)
        return 2

    if args.format == "CSV":
        write_csv(records, args.output)
    else:
        write_fasta(
            records,
            args.output,
            annotate=args.fasta_annotate,
            include_names=args.fasta_name,
        )

    print(f"Wrote {len(records)} records to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

