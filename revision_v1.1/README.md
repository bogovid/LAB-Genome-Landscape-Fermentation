# FEMS revision materials (v1.1.0)

This directory records the validation and genome-wide analyses added for the FEMS-ready revision of the LAB Genome Landscape manuscript.

## Included in GitHub

- `validation/FM02_FM06_sequence_validation.csv`: BLASTP/TBLASTN validation summary for critical FM02/FM06 non-detections.
- `validation/gene_call_corrections.csv`: validated gene-call corrections; currently one correction, LABG_0506 GC025/alsD from 0 to 1.
- `genomewide/genomewide_summary.txt`: pairwise Mash/FastANI summary used in the manuscript.

## Full versioned data archive

The complete v1.1.0 revision data package is intended for the corresponding Zenodo version and contains, in addition to the files above:

- corrected 507-genome fermentation marker matrix;
- recomputed assembly contiguity metrics for all genomes;
- collapsed 28-marker profiles;
- recomputed manuscript-supporting summary tables;
- gzip-compressed unique-pair Mash/marker and FastANI/marker tables.

Raw reciprocal all-vs-all Mash and FastANI outputs are not duplicated because the unique-pair tables preserve every value used in the manuscript analyses while removing self-comparisons and reciprocal duplicates.

## Interpretation

Marker presence denotes detected genomic coding potential within the operational marker framework. It is not evidence of expression, enzyme activity, metabolic flux, or fermentation phenotype.
