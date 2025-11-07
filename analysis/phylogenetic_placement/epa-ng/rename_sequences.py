#!/usr/bin/env python3
"""
Rename FASTA sequences to simple IDs (ref001, ref002, etc.) for EPA-ng/RAxML-ng compatibility.
Creates a mapper and saves both the renamed FASTA and JSON mapper.
"""

import json
import sys
from pathlib import Path
from Bio import SeqIO


def create_simple_mapper(fasta_path, output_fasta, mapper_json, prefix="ref"):
    """
    Create simple sequence IDs and save mapper.

    Parameters:
    -----------
    fasta_path : str or Path
        Input FASTA file path
    output_fasta : str or Path
        Output FASTA file with renamed sequences
    mapper_json : str or Path
        Output JSON file with ID mapping
    prefix : str
        Prefix for simple IDs (default: "ref")
    """
    fasta_path = Path(fasta_path)
    output_fasta = Path(output_fasta)
    mapper_json = Path(mapper_json)

    # Create output directory if needed
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    mapper_json.parent.mkdir(parents=True, exist_ok=True)

    # Read sequences and create mapper
    mapper = {}
    renamed_records = []

    print(f"Reading sequences from: {fasta_path}")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Found {len(records)} sequences")

    for idx, record in enumerate(records, start=1):
        # Create simple ID
        simple_id = f"{prefix}{idx:03d}"

        # Store mapping (both directions for convenience)
        mapper[simple_id] = {"original_id": record.id, "original_description": record.description}

        # Create renamed record
        record.id = simple_id
        record.name = simple_id
        record.description = f"{simple_id} (original: {record.description})"
        renamed_records.append(record)

    # Write renamed FASTA
    print(f"Writing renamed FASTA to: {output_fasta}")
    SeqIO.write(renamed_records, output_fasta, "fasta")

    # Also create reverse mapper for easy lookup
    reverse_mapper = {v["original_id"]: k for k, v in mapper.items()}

    full_mapper = {
        "simple_to_original": mapper,
        "original_to_simple": reverse_mapper,
        "metadata": {
            "total_sequences": len(records),
            "prefix": prefix,
            "input_file": str(fasta_path),
            "output_file": str(output_fasta),
        },
    }

    # Write mapper JSON
    print(f"Writing mapper to: {mapper_json}")
    with open(mapper_json, "w") as f:
        json.dump(full_mapper, f, indent=2)

    print(f"\nSuccess!")
    print(f"  - Renamed {len(records)} sequences")
    print(f"  - Output FASTA: {output_fasta}")
    print(f"  - Mapper JSON: {mapper_json}")

    return full_mapper


def main():
    """Command-line interface"""
    if len(sys.argv) < 3:
        print("Usage: python rename_sequences.py <input.fasta> <output.fasta> [mapper.json] [prefix]")
        print("\nExample:")
        print("  python rename_sequences.py reference.fasta reference_renamed.fasta mapper.json ref")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    mapper_json = sys.argv[3] if len(sys.argv) > 3 else "sequence_mapper.json"
    prefix = sys.argv[4] if len(sys.argv) > 4 else "ref"

    create_simple_mapper(input_fasta, output_fasta, mapper_json, prefix)


if __name__ == "__main__":
    main()
