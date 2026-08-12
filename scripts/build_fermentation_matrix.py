#!/usr/bin/env python3

import csv
import re
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]

CATALOG = PROJECT / "metadata" / "processed" / "gene_catalog_analysis_ready.csv"
METADATA = PROJECT / "metadata" / "processed" / "genome_metadata.csv"
ANNOT_ROOT = PROJECT / "analysis" / "prokka_annotations"

OUTDIR = PROJECT / "analysis" / "functional_modules" / "fermentation"
OUTDIR.mkdir(parents=True, exist_ok=True)

HITS_FILE = OUTDIR / "fermentation_marker_hits.tsv"
GENE_MATRIX = OUTDIR / "fermentation_gene_matrix.csv"
MODULE_SUMMARY = OUTDIR / "fermentation_module_summary.csv"


def normalize_gene(gene):
    """
    Prokka may append suffixes such as _1, _2 to duplicated genes.
    Example: ldh_1 -> ldh
    """
    gene = (gene or "").strip()
    return re.sub(r"_\d+$", "", gene)


def split_ec(ec_text):
    """
    Convert EC strings into a comparable set.
    Handles forms such as:
    1.1.1.1
    1.2.1.10/1.1.1.1
    1.2.1.10;1.1.1.1
    """
    if not ec_text:
        return set()

    return {
        x.strip()
        for x in re.split(r"[/;, ]+", ec_text.strip())
        if x.strip()
    }


def load_catalog():
    markers = []

    with CATALOG.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            markers.append({
                "Gene_ID": row["Gene_ID"].strip(),
                "Module_ID": row["Module_ID"].strip(),
                "Gene_symbol": row["Gene_symbol"].strip(),
                "Gene_symbol_norm": normalize_gene(row["Gene_symbol"]),
                "EC_numbers": split_ec(row["EC_number"]),
                "Enzyme_name": row["Enzyme_name"].strip(),
            })

    return markers


def load_metadata():
    with METADATA.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def marker_match(feature, marker):
    gene_raw = feature.get("gene", "")
    gene_norm = normalize_gene(gene_raw)

    feature_ecs = split_ec(feature.get("EC_number", ""))

    gene_match = (
        bool(marker["Gene_symbol_norm"])
        and gene_norm.lower() == marker["Gene_symbol_norm"].lower()
    )

    ec_match = bool(
        marker["EC_numbers"]
        and feature_ecs
        and marker["EC_numbers"].intersection(feature_ecs)
    )

    if gene_match and ec_match:
        evidence = "gene+EC"
    elif gene_match:
        evidence = "gene"
    elif ec_match:
        evidence = "EC"
    else:
        return None

    return evidence


def main():
    markers = load_catalog()
    genomes = load_metadata()

    print(f"Genoma u metadata tabeli: {len(genomes)}")
    print(f"Markera u analitičkom katalogu: {len(markers)}")
    print()

    marker_ids = [m["Gene_ID"] for m in markers]

    marker_to_module = {
        m["Gene_ID"]: m["Module_ID"]
        for m in markers
    }

    module_markers = defaultdict(list)
    for marker in markers:
        module_markers[marker["Module_ID"]].append(marker["Gene_ID"])

    presence = {}
    all_hits = []

    missing_tsv = []

    for index, genome in enumerate(genomes, start=1):
        genome_id = genome["Genome_ID"]
        species = genome["Species_current"]
        accession = genome["Accession"]

        tsv = ANNOT_ROOT / genome_id / f"{genome_id}.tsv"

        genome_presence = {marker_id: 0 for marker_id in marker_ids}

        if not tsv.is_file():
            missing_tsv.append(genome_id)
            presence[genome_id] = genome_presence
            continue

        with tsv.open(encoding="utf-8", newline="") as handle:
            features = list(csv.DictReader(handle, delimiter="\t"))

        for feature in features:
            if feature.get("ftype", "") != "CDS":
                continue

            for marker in markers:
                evidence = marker_match(feature, marker)

                if evidence is None:
                    continue

                genome_presence[marker["Gene_ID"]] = 1

                all_hits.append({
                    "Genome_ID": genome_id,
                    "Species": species,
                    "Accession": accession,
                    "Gene_ID": marker["Gene_ID"],
                    "Module_ID": marker["Module_ID"],
                    "Catalog_gene_symbol": marker["Gene_symbol"],
                    "locus_tag": feature.get("locus_tag", ""),
                    "Prokka_gene": feature.get("gene", ""),
                    "EC_number": feature.get("EC_number", ""),
                    "product": feature.get("product", ""),
                    "Evidence": evidence,
                })

        presence[genome_id] = genome_presence

        if index % 25 == 0 or index == len(genomes):
            print(f"Obrađeno: {index}/{len(genomes)}")

    # ---------------------------------------------------------
    # 1. Detailed hits
    # ---------------------------------------------------------

    hit_fields = [
        "Genome_ID",
        "Species",
        "Accession",
        "Gene_ID",
        "Module_ID",
        "Catalog_gene_symbol",
        "locus_tag",
        "Prokka_gene",
        "EC_number",
        "product",
        "Evidence",
    ]

    with HITS_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=hit_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(all_hits)

    # ---------------------------------------------------------
    # 2. Gene presence / absence matrix
    # ---------------------------------------------------------

    gene_fields = [
        "Genome_ID",
        "Species",
        "Accession",
        "Study_group",
    ] + marker_ids

    with GENE_MATRIX.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_fields)
        writer.writeheader()

        for genome in genomes:
            genome_id = genome["Genome_ID"]

            row = {
                "Genome_ID": genome_id,
                "Species": genome["Species_current"],
                "Accession": genome["Accession"],
                "Study_group": genome["Study_group"],
            }

            row.update(presence[genome_id])
            writer.writerow(row)

    # ---------------------------------------------------------
    # 3. Module coverage summary
    # ---------------------------------------------------------

    module_ids = sorted(module_markers)

    module_fields = [
        "Genome_ID",
        "Species",
        "Accession",
        "Study_group",
    ]

    for module_id in module_ids:
        module_fields.extend([
            f"{module_id}_present",
            f"{module_id}_total",
            f"{module_id}_coverage",
        ])

    with MODULE_SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=module_fields)
        writer.writeheader()

        for genome in genomes:
            genome_id = genome["Genome_ID"]

            row = {
                "Genome_ID": genome_id,
                "Species": genome["Species_current"],
                "Accession": genome["Accession"],
                "Study_group": genome["Study_group"],
            }

            for module_id in module_ids:
                ids = module_markers[module_id]
                present = sum(presence[genome_id][gene_id] for gene_id in ids)
                total = len(ids)

                row[f"{module_id}_present"] = present
                row[f"{module_id}_total"] = total
                row[f"{module_id}_coverage"] = round(
                    present / total,
                    3,
                )

            writer.writerow(row)

    print()
    print(f"Ukupno pogodaka: {len(all_hits)}")
    print(f"TSV fajlova koji nedostaju: {len(missing_tsv)}")
    print()
    print(f"Hits: {HITS_FILE}")
    print(f"Gene matrix: {GENE_MATRIX}")
    print(f"Module summary: {MODULE_SUMMARY}")


if __name__ == "__main__":
    main()
