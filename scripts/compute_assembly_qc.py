#!/usr/bin/env python3
"""Compute assembly contiguity metrics from genome FASTA files.

Usage:
    python scripts/compute_assembly_qc.py --genomes-dir PATH --output genome_assembly_qc.csv

The directory may contain .fna files directly or recursively.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_lengths(path: Path) -> list[int]:
    lengths = []
    current = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
        if current:
            lengths.append(current)
    return lengths


def n50(lengths: list[int]) -> int:
    total = sum(lengths)
    half = total / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def genome_id_from_path(path: Path) -> str:
    stem = path.stem
    if stem.startswith("LABG_"):
        return stem
    for part in path.parts:
        if part.startswith("LABG_"):
            return part
    return stem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    files = sorted(args.genomes_dir.rglob("*.fna"))
    rows = []
    for path in files:
        lengths = fasta_lengths(path)
        rows.append(
            {
                "Genome_ID": genome_id_from_path(path),
                "Contigs_recomputed": len(lengths),
                "Genome_size_recomputed_bp": sum(lengths),
                "N50_recomputed_bp": n50(lengths),
                "Longest_contig_bp": max(lengths) if lengths else 0,
                "QC_status": "OK" if lengths else "EMPTY_FASTA",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [
            "Genome_ID", "Contigs_recomputed", "Genome_size_recomputed_bp",
            "N50_recomputed_bp", "Longest_contig_bp", "QC_status"
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Genomes: {len(rows)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
