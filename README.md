# LAB Genome Landscape: Fermentation-associated genomic architecture in lactic acid bacteria

This repository contains the code and compact reproducibility materials supporting a comparative genomic analysis of fermentation-associated genomic markers in a curated collection of 507 bacterial genomes.

## Dataset

- 507 total genomes
- 486 core lactic acid bacteria (LAB)
- 10 fermentative non-LAB genomes
- 8 external comparators
- 3 Enterococcus genomes
- 462 Lactiplantibacillus genomes

The original RefSeq FASTA assemblies are not redistributed. Public assembly accession numbers are provided in `metadata/refseq_accession_manifest.csv`.

## Functional framework

Seven fermentation-associated modules (FM01-FM07) were analysed using 28 unique biological markers derived from 31 validated catalogue entries.

- FM01: EMP glycolysis
- FM02: Lactate formation
- FM03: Phosphoketolase / pentose-phosphate-associated metabolism
- FM04: Acetate branch
- FM05: Ethanol branch
- FM06: Alternative pyruvate fates
- FM07: NAD/redox balance

Marker presence is interpreted as genomic coding potential within this operational framework, not as direct evidence of gene expression, metabolic flux, enzyme activity, or fermentation phenotype.

## v1.1.0 FEMS revision

Version 1.1.0 adds validation and genome-wide analyses used in the FEMS-ready revision of the manuscript.

### Sequence validation

Critical apparent absences in FM02 and FM06 were independently checked with BLASTP and TBLASTN. One Prokka-based false-negative was identified and corrected:

- LABG_0506 (Weissella confusa), GC025 / alsD: 0 -> 1.

For the three Lactiplantibacillus genomes lacking pfl in the annotation-based matrix, no BLASTP or TBLASTN homolog was detected, supporting the retained non-detection.

### Assembly quality

Assembly contiguity metrics were recomputed from all 507 genome FASTA files. The three Lactiplantibacillus pfl-negative genomes are highly contiguous (3-4 contigs; N50 approximately 2.45-3.04 Mb), arguing against fragmentation as a simple explanation for the observed non-detection.

### Genome-wide comparison

For the 486 core-LAB genomes, 117,855 unique Mash genome pairs were compared with collapsed 28-marker profiles:

- 106,047 pairs (89.98%) had identical 28-marker profiles.
- all 95,035 pairs with Mash distance <=0.10 had identical profiles.
- Spearman rho between Mash distance and marker discordance = 0.514690.

For the 462 Lactiplantibacillus genomes, 106,491 unique FastANI pairs were analysed:

- ANI range: 79.4657-100.0000%.
- 105,114 pairs (98.71%) had identical 28-marker profiles.
- the remaining 1,377 pairs differed at exactly one marker.

These results support the conclusion that the targeted 28-marker repertoire remains highly conserved across substantial whole-genome divergence.

## Reproducibility layout

The repository keeps code and compact validation summaries. FEMS-revision material is under `revision_v1.1/`, including:

- `validation/`: sequence-level validation and gene-call correction log;
- `genomewide/`: genome-wide comparison summary;
- `README.md`: description of the complete versioned archive.

The complete v1.1.0 data snapshot is archived separately on Zenodo and includes the corrected 507-genome marker matrix, assembly-QC table, collapsed 28-marker profiles, recomputed summary tables, and gzip-compressed unique-pair Mash/FastANI tables. Raw reciprocal all-vs-all Mash and FastANI outputs are omitted because the unique-pair tables preserve all values used in the manuscript while removing self-comparisons and reciprocal duplicates.

## Archived releases

- Software/repository snapshot v1.1.0: DOI `10.5281/zenodo.22718746`
- Supporting data snapshot: DOI `10.5281/zenodo.22718859`
- GitHub release: `v1.1.0`

## Reproduction

Run `scripts/analyze_fermentation_landscape.py` to regenerate the principal downstream tables.
Run `scripts/plot_genus_module_heatmap.py` to regenerate Figure 1.
Run `scripts/compute_assembly_qc.py` to recompute assembly contiguity metrics from FASTA files.
Run `scripts/summarize_genomewide_marker_conservation.py` on the versioned unique-pair tables to reproduce the principal pairwise summaries.

Software versions are recorded in `reproducibility/software_versions.txt`. The v1.1 validation used BLAST+ 2.16.0+, Mash 2.3, and FastANI 1.33.

## Version

Release: v1.1.0
Dataset snapshot: LAB Genome Landscape v1.1
