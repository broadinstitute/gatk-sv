# 02. Evidence assessment

This module provides workflows for assessing read-depth (RD), paired-end (PE),
split read (SR), and B-allele frequency (BAF) evidence supporting structural
variant calls.

## Required matrices
To assess the evidence for potential SVs, there are three matrices required as input:
* `Matrics.binCov.bed.gz` and `Matrics.binCov.median`
* `Matrics.pe.sorted.txt.gz`
* `Matrics.sr.sorted.txt.gz`

All matrices should have been properly bgzipped and tabix indexed. Refer to each subdirectory for details of how these matrics are created and applied to the analysis.

For further details for processing each evidence, please refer to the readme under each sub-directory.
