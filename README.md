# LAB Genome Landscape: Fermentation-associated genomic architecture in lactic acid bacteria

This repository contains the reproducibility package supporting a comparative genomic analysis of fermentation-associated functional modules in a curated collection of 507 bacterial genomes.

## Dataset

- 507 total genomes
- 486 core lactic acid bacteria (LAB)
- 10 fermentative non-LAB genomes
- 8 external comparators
- 3 Enterococcus genomes
- 462 Lactiplantibacillus genomes

The original RefSeq FASTA assemblies are not redistributed. Public assembly accession numbers are provided in metadata/refseq_accession_manifest.csv.

## Functional framework

Seven fermentation-associated modules (FM01-FM07) were analysed using 28 unique biological markers derived from 31 validated catalogue entries.

- FM01: EMP glycolysis
- FM02: Lactate formation
- FM03: Phosphoketolase / pentose-phosphate-associated metabolism
- FM04: Acetate branch
- FM05: Ethanol branch
- FM06: Alternative pyruvate fates
- FM07: NAD/redox balance

## Main reproducible findings

- 27/28 markers were detected in all 462 Lactiplantibacillus genomes.
- pfl was detected in 459/462 genomes (99.35%).
- All 28 markers were detected in all genomes of L. plantarum (n=411), L. pentosus (n=21), L. argentoratensis (n=16), and L. paraplantarum (n=9).
- FM04 showed complete mean marker representation across sampled multi-genome LAB genera.
- FM02 showed the largest inter-genus prevalence range.
- FM06 showed the highest coefficient of variation.

These results describe genomic capacity within the standardized marker framework and do not constitute direct measurements of metabolic flux or fermentation phenotype.

## Reproduction

Run scripts/analyze_fermentation_landscape.py to regenerate the principal downstream tables.
Run scripts/plot_genus_module_heatmap.py to regenerate Figure 1.

Software versions are recorded in reproducibility/software_versions.txt.

## Version

Release: v1.0
Dataset snapshot: LAB Genome Landscape v1.0
