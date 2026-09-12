# LAB Genome Landscape Fermentation v1.1.0

This release supports the FEMS-ready revision of the manuscript *Conserved fermentation-associated genomic architecture and lineage-associated variation in lactic acid bacteria*.

## What changed since v1.0.0

- Added independent BLASTP/TBLASTN validation of critical FM02/FM06 apparent absences.
- Corrected one annotation false-negative: LABG_0506 (*Weissella confusa*), GC025/`alsD`, from 0 to 1.
- Added assembly-contiguity QC based on recomputed contig counts and N50 values for all 507 genomes.
- Added genome-wide Mash analysis for 486 core-LAB genomes.
- Added all-vs-all FastANI analysis for 462 *Lactiplantibacillus* genomes.
- Added collapsed 28-marker pairwise conservation analysis.
- Updated interpretation to distinguish genomic marker presence from expression, metabolic flux, enzyme activity, or fermentation phenotype.

## Key v1.1.0 results

- 117,855 unique core-LAB Mash pairs were analysed; 106,047 (89.98%) had identical 28-marker profiles.
- All 95,035 pairs with Mash distance <=0.10 had identical 28-marker profiles.
- Spearman rho between Mash distance and marker discordance was 0.514690.
- 106,491 unique *Lactiplantibacillus* FastANI pairs were analysed across ANI values of 79.4657-100.0000%.
- 105,114 pairs (98.71%) had identical 28-marker profiles; the remaining 1,377 pairs differed at exactly one marker.
- The three independently validated `pfl`-negative *Lactiplantibacillus* genomes retained their non-detection and were not explained by poor assembly contiguity.

## Data archiving

The GitHub repository contains code and compact validation summaries. The complete versioned data snapshot is archived on Zenodo and includes the corrected marker matrix, assembly-QC table, collapsed 28-marker profiles, recomputed summary tables, and gzip-compressed unique-pair Mash/FastANI tables.

The underlying NCBI RefSeq genome assemblies are not redistributed.
