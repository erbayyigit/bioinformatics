#!/usr/bin/env python3
"""Mutate a coding DNA sequence by amino-acid position.

Usage example:
  python3 scripts/02_mutate_coding_sequence.py -m ARG2ALA
  python3 scripts/02_mutate_coding_sequence.py sample_data/02_mutate_coding_sequence_input.fasta -m E2A
  python3 scripts/02_mutate_coding_sequence.py -m K3A,E2A
  python3 scripts/02_mutate_coding_sequence.py -m [K3A,E2A]

The script accepts one or more `-m`/`--mutation` flags. Each mutation has the
form `ARG2ALA` or `R2A` (three-letter or one-letter AA codes). Simple mutations
separated by commas are written as separate records. Bracketed groups like
`[K3A,E2A]` are treated as a compound mutation applied together. Positions are
1-based amino-acid positions. The script verifies the original amino-acid at
the requested position and writes a new FASTA with the mutated DNA; the
modified codon is written in lower-case letters. Output filename includes the
record id and mutations.
"""

import re
import sys
from pathlib import Path
import argparse
from typing import List, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Data import CodonTable
except Exception:
    SeqIO = None
    Seq = None
    SeqRecord = None
    CodonTable = None


# Default output directory: sibling `output/` next to `src/` inside the project
ROOT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ROOT_DIR / "output"
DEFAULT_INPUT = ROOT_DIR / "sample_data" / "02_mutate_coding_sequence_input.fasta"

STANDARD_DNA_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


ecoli_codon_table = {
    '*': {'UAA': 0.64, 'UAG': 0.07, 'UGA': 0.29},
    'A': {'GCA': 0.21, 'GCC': 0.27, 'GCG': 0.36, 'GCU': 0.16},
    'C': {'UGC': 0.56, 'UGU': 0.44},
    'D': {'GAC': 0.37, 'GAU': 0.63},
    'E': {'GAA': 0.69, 'GAG': 0.31},
    'F': {'UUC': 0.43, 'UUU': 0.57},
    'G': {'GGA': 0.11, 'GGC': 0.41, 'GGG': 0.15, 'GGU': 0.34},
    'H': {'CAC': 0.43, 'CAU': 0.57},
    'I': {'AUA': 0.07, 'AUC': 0.42, 'AUU': 0.51},
    'K': {'AAA': 0.76, 'AAG': 0.24},
    'L': {'CUA': 0.04, 'CUC': 0.1, 'CUG': 0.5, 'CUU': 0.1, 'UUA': 0.13, 'UUG': 0.13},
    'M': {'AUG': 1.0},
    'N': {'AAC': 0.55, 'AAU': 0.45},
    'P': {'CCA': 0.19, 'CCC': 0.12, 'CCG': 0.53, 'CCU': 0.16},
    'Q': {'CAA': 0.35, 'CAG': 0.65},
    'R': {'AGA': 0.04, 'AGG': 0.02, 'CGA': 0.06, 'CGC': 0.4, 'CGG': 0.1, 'CGU': 0.38},
    'S': {'AGC': 0.28, 'AGU': 0.15, 'UCA': 0.12, 'UCC': 0.15, 'UCG': 0.15, 'UCU': 0.15},
    'T': {'ACA': 0.13, 'ACC': 0.44, 'ACG': 0.27, 'ACU': 0.16},
    'V': {'GUA': 0.15, 'GUC': 0.22, 'GUG': 0.37, 'GUU': 0.26},
    'W': {'UGG': 1.0},
    'Y': {'UAC': 0.43, 'UAU': 0.57}
}


# Minimal three-letter to one-letter mapping for standard amino acids
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def normalize_aa(code: str) -> str:
    """Normalize an amino-acid code (one- or three-letter) to one-letter uppercase."""
    code = code.strip().upper()
    if len(code) == 1:
        return code
    if len(code) == 3:
        return THREE_TO_ONE.get(code, "?")
    return "?"


def parse_mutation(token: str) -> Tuple[str, int, str]:
    """Parse a mutation token like ARG2ALA or R2A into (orig_aa, position, new_aa)."""
    m = re.match(r"^([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})$", token.strip())
    if not m:
        raise ValueError(f"Invalid mutation format: {token}")
    a1, pos, a2 = m.groups()
    o = normalize_aa(a1)
    n = normalize_aa(a2)
    if o == "?" or n == "?":
        raise ValueError(f"Unknown amino-acid code in mutation: {token}")
    return o, int(pos), n


def parse_mutation_groups(tokens: List[str]) -> List[List[Tuple[str, int, str]]]:
    """Parse mutation groups from one or more -m values.

    - Simple mutations separated by commas are treated as independent records.
    - Bracketed groups like [K3A,E2A] are treated as a compound mutation.
    """
    groups = []
    for token in tokens:
        token = token.strip()
        if not token:
            continue

        if token.startswith("[") and token.endswith("]"):
            inner = token[1:-1].strip()
            if not inner:
                raise ValueError(f"Compound mutation group is empty: {token}")
            parts = [part.strip() for part in inner.split(",") if part.strip()]
            if not parts:
                raise ValueError(f"Compound mutation group is empty: {token}")
            groups.append([parse_mutation(part) for part in parts])
        else:
            parts = [part.strip() for part in token.split(",") if part.strip()]
            for part in parts:
                groups.append([parse_mutation(part)])
    return groups


def mutation_label(mutation: Tuple[str, int, str]) -> str:
    orig_aa, pos, new_aa = mutation
    return f"{orig_aa}{pos}{new_aa}"


def translate_codon(codon: str) -> str:
    return STANDARD_DNA_CODON_TABLE.get(codon.upper(), None)


def translate_dna(dna_seq: str) -> str:
    if len(dna_seq) % 3 != 0:
        raise ValueError("Input DNA length is not a multiple of 3")
    amino_acids = []
    for start in range(0, len(dna_seq), 3):
        codon = dna_seq[start:start + 3].upper()
        aa = translate_codon(codon)
        if aa is None:
            raise ValueError(f"Cannot translate codon '{codon}'")
        amino_acids.append(aa)
    return "".join(amino_acids)


def validate_mutation_group(dna_seq: str, mutations: List[Tuple[str, int, str]]) -> None:
    """Validate a compound mutation group against the original DNA sequence."""
    seq = dna_seq.upper()
    if len(seq) % 3 != 0:
        raise ValueError("Input DNA length is not a multiple of 3")

    seen_positions = set()

    for orig_aa, pos, _ in mutations:
        if pos in seen_positions:
            raise ValueError(f"Duplicate mutation position specified within group: {pos}")
        seen_positions.add(pos)

        codon_start = (pos - 1) * 3
        codon = seq[codon_start:codon_start + 3]
        if len(codon) != 3:
            raise ValueError(f"Mutation position {pos} is out of range for sequence length")

        aa = translate_codon(codon)
        if aa is None:
            raise ValueError(f"Cannot translate codon '{codon}' at position {pos}")
        if aa != orig_aa:
            raise ValueError(f"Mismatch at position {pos}: expected {orig_aa}, found {aa}")


def build_reverse_codon_map() -> dict:
    rev = {}
    for codon, aa in STANDARD_DNA_CODON_TABLE.items():
        rev.setdefault(aa, []).append(codon.upper())
    rev["*"] = [c.upper() for c in STOP_CODONS]
    return rev


def apply_mutation(dna_seq: str, mutation: Tuple[str, int, str], rev_map: dict) -> str:
    """Apply a single mutation to the DNA sequence and return mutated DNA string.

    The mutated codon is converted to lower-case; other sequence remains upper-case.
    Position is 1-based amino-acid position.
    """
    orig_aa, pos, new_aa = mutation
    seq = dna_seq.upper()
    if len(seq) % 3 != 0:
        raise ValueError("Input DNA length is not a multiple of 3")
    codon_start = (pos - 1) * 3
    codon = seq[codon_start:codon_start + 3]
    if len(codon) != 3:
        raise ValueError(f"Position {pos} is out of range for sequence length")
    # translate codon to AA
    aa = translate_codon(codon)
    if aa is None:
        raise ValueError(f"Cannot translate codon '{codon}' at position {pos}")
    if aa != orig_aa:
        raise ValueError(f"Mismatch at position {pos}: expected {orig_aa}, found {aa}")

    # choose replacement codon for new_aa using E. coli frequency data
    candidates = rev_map.get(new_aa, [])
    if not candidates:
        raise ValueError(f"No codons found for amino-acid '{new_aa}'")

    # sort candidate codons by E. coli frequency descending
    def codon_frequency(candidate: str) -> float:
        rna_codon = candidate.replace("T", "U")
        return ecoli_codon_table.get(new_aa, {}).get(rna_codon, 0.0)

    sorted_candidates = sorted(candidates, key=codon_frequency, reverse=True)
    replacement = sorted_candidates[0]

    # if the highest-frequency codon is the same as the existing codon, prefer the next highest-frequency codon if available
    if replacement.upper() == codon.upper() and len(sorted_candidates) > 1:
        replacement = sorted_candidates[1]

    # build mutated sequence with replacement codon in lower-case
    mutated = seq[:codon_start] + replacement.lower() + seq[codon_start + 3:]
    return mutated


def apply_mutation_group(dna_seq: str, mutations: List[Tuple[str, int, str]], rev_map: dict) -> str:
    """Apply a list of mutations to the DNA sequence as a group."""
    mutated = dna_seq
    for mutation in mutations:
        mutated = apply_mutation(mutated, mutation, rev_map)
    return mutated


def read_fasta_records(file_path: str):
    if SeqIO is not None:
        try:
            return list(SeqIO.parse(file_path, "fasta"))
        except Exception:
            pass

    class SimpleRecord:
        def __init__(self, record_id: str, sequence: str):
            self.id = record_id
            self.seq = sequence

    records = []
    current_id = None
    current_seq = []
    with open(file_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append(SimpleRecord(current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append(SimpleRecord(current_id, "".join(current_seq)))
    return records


def write_fasta_records(file_path: Path, records: List[Tuple[str, str, str]]) -> None:
    with open(file_path, "w", encoding="utf-8") as handle:
        for record_id, description, sequence in records:
            handle.write(f">{record_id}")
            if description:
                handle.write(f" {description}")
            handle.write("\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx:idx + 80] + "\n")


def main(argv: List[str] = None) -> int:
    parser = argparse.ArgumentParser(description="Apply amino-acid mutations to a coding FASTA record")
    parser.add_argument(
        "input_fasta",
        nargs="?",
        default=str(DEFAULT_INPUT),
        help="Input FASTA file (single record). Defaults to sample_data/02_mutate_coding_sequence_input.fasta",
    )
    parser.add_argument(
        "-m",
        "--mutation",
        action="append",
        required=True,
        help="Mutation token, e.g. ARG2ALA or R2A. Separate simple mutations with commas or repeated -m flags; use [K3A,E2A] for a compound mutation.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args(argv)

    input_arg = args.input_fasta
    input_path = Path(input_arg)
    if not input_path.is_absolute():
        cwd_candidate = Path.cwd() / input_path
        if cwd_candidate.exists():
            input_file = str(cwd_candidate.resolve())
        else:
            input_file = str((ROOT_DIR / input_path).resolve())
    else:
        input_file = str(input_path)
    args.input_fasta = input_file

    if not Path(input_file).exists():
        print(f"Error: input FASTA file not found: {input_file}", file=sys.stderr)
        return 2

    records = read_fasta_records(input_file)
    if len(records) != 1:
        print("Error: input FASTA must contain exactly one record", file=sys.stderr)
        return 2
    rec = records[0]
    dna = str(rec.seq).upper()

    # basic validation: start with ATG and end with stop codon
    if not dna.startswith("ATG"):
        print("Error: sequence does not start with ATG", file=sys.stderr)
        return 3
    if dna[-3:].upper() not in STOP_CODONS:
        print("Error: sequence does not end with a stop codon", file=sys.stderr)
        return 4

    # build rev map
    rev_map = build_reverse_codon_map()

    # parse mutation groups and validate against the original sequence
    parsed_groups = parse_mutation_groups(args.mutation)
    try:
        for group in parsed_groups:
            validate_mutation_group(dna, group)
    except ValueError as exc:
        print(f"Error validating mutations: {exc}", file=sys.stderr)
        return 5

    # apply each group separately to the original sequence
    mutated_groups = []
    try:
        for group in parsed_groups:
            mutated_groups.append(apply_mutation_group(dna, group, rev_map))
    except ValueError as exc:
        print(f"Error applying mutation: {exc}", file=sys.stderr)
        return 5

    # Check that each mutated sequence does not introduce a premature stop codon
    try:
        orig_protein = translate_dna(dna)
        for mutated_seq in mutated_groups:
            mut_protein = translate_dna(mutated_seq)
            orig_first_stop = orig_protein.find("*")
            mut_first_stop = mut_protein.find("*")
            if orig_first_stop == -1:
                orig_first_stop = len(orig_protein)
            if mut_first_stop == -1:
                mut_first_stop = len(mut_protein)
            if mut_first_stop < orig_first_stop:
                print(
                    f"Error: mutation introduces a premature stop codon at protein position {mut_first_stop + 1}",
                    file=sys.stderr,
                )
                return 6
    except Exception as exc:
        print(f"Error checking mutated protein sequence: {exc}", file=sys.stderr)
        return 6

    # prepare output
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    mut_suffix = "_".join(
        mutation_label(mut).lower() for group in parsed_groups for mut in group
    )
    out_name = f"{rec.id}_{mut_suffix}.fasta"
    out_path = out_dir / out_name

    # write DNA FASTA with original DNA first, then one mutated DNA record per parsed group
    dna_records = [(f"{rec.id}_dna", "Original DNA coding sequence", dna)]
    mutated_dna_seqs = list(mutated_groups)
    for group in parsed_groups:
        label = "_".join(mutation_label(mut) for mut in group)
        mutated_seq = mutated_dna_seqs.pop(0)
        dna_records.append(
            (f"{rec.id}_{label}", f"Mutated DNA coding sequence; Mutations: {', '.join(mutation_label(mut) for mut in group)}", mutated_seq)
        )
    write_fasta_records(out_path, dna_records)
    print(f"Wrote DNA FASTA: {out_path}")

    # write protein FASTA with original protein first, then one mutated protein record per parsed group
    orig_protein = translate_dna(dna)
    prot_out_name = f"{rec.id}_{mut_suffix}_prot.fasta"
    prot_out_path = out_dir / prot_out_name
    prot_records = [(f"{rec.id}_protein", "Original protein sequence", orig_protein)]
    mutated_dna_seqs = list(mutated_groups)
    for group in parsed_groups:
        label = "_".join(mutation_label(mut) for mut in group)
        mutated_protein = translate_dna(mutated_dna_seqs.pop(0))
        prot_records.append(
            (f"{rec.id}_{label}_prot", f"Mutated protein sequence; Mutations: {', '.join(mutation_label(mut) for mut in group)}", mutated_protein)
        )
    write_fasta_records(prot_out_path, prot_records)
    print(f"Wrote protein FASTA: {prot_out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
