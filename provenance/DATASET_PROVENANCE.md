# Dataset Provenance

## Project

LAB Genome Landscape v1.0

Working title:

Comparative Genomic Landscape of a Curated Reference Collection of Lactic Acid Bacteria

---

## Primary data source

National Center for Biotechnology Information (NCBI)

Reference genome assemblies downloaded from the RefSeq collection.

---

## Dataset construction

The genome collection was assembled through a reproducible curation workflow.

Main steps included:

- retrieval of RefSeq genome assemblies
- exclusion of duplicate genome records
- preference for RefSeq assemblies over duplicate GenBank entries
- organization into a local reference genome library
- standardized metadata generation
- genome path validation
- creation of a frozen project snapshot (v1.0)

---

## Frozen dataset

Snapshot version:

LAB Genome Landscape v1.0

Frozen before any comparative genomic analyses.

Dataset size:

507 genome assemblies

All genome files successfully validated.

---

## Metadata

The primary metadata table is

metadata/processed/genome_metadata.csv

This table serves as the single source of truth for all downstream analyses.

---

## FASTA integrity

A complete genome manifest and SHA-256 checksums were generated for every FASTA file.

Files:

manifests/fasta_manifest_v1.0.tsv

checksums/fasta_sha256_v1.0.tsv

---

## Reproducibility

All comparative genomic analyses performed in this project must explicitly reference:

LAB Genome Landscape v1.0

Later updates of the genome library must be released as new snapshot versions (e.g. v1.1, v2.0) and must not overwrite this frozen dataset.

Created by: Bogovid Živković
Project: LAB Genome Landscape
