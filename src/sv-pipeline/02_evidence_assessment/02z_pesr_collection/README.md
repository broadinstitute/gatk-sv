# pesrTest
This repository contains the workflow that combines and evaluates pair-end and split reads support for all SV calls on a per-batch basis.

## Process
To evaluate the pair-end and split read evidence, the following process should be applied to each sample:
1. count discordant read pairs and split reads:
```
python scripts/count_disc.py -s sample_name -r chrom sample_name.bam |sort -k1,1V -k2,2n -k5,5n|bgzip -c > disc_counts/chrom/sample_name.txt.gz
python scripts/count_disc.py --tloc -s sample_name -r chrom sample_name.bam|sort -k1,1V -k2,2n -k5,5n|bgzip -c > tloc_counts/chrom/sample_name.txt.gz
python scripts/count_splits.py -r chrom sample_name.bam stdout | bgzip -c > split_counts/chrom/sample_name.txt.gz
```
2. Filter each of the discordant information with whitelist excluding the genomic regions of decent mappability, as the following example.  `ref/b37.lumpy.include.4-13.bed.gz` is provided as an example of whitelist for hs37, while the proper list for different reference should have been prepared accordingly.
```
tabix -s1 -b2 -e2 disc_counts/chrom/sample_name.txt.gz
tabix -R whitelist disc_counts/chrom/sample_name.txt.gz | bgzip -c > disc_counts_filtered/chrom/sample_name.txt.gz
tabix -s1 -b2 -e2 disc_counts_filtered/chrom/sample_name.txt.gz
```
3. merge disc_counts and tloc_counts:
```
sort -m -k1,1V -k2,2n -k4,4V -k5,5n <(bgzip -d -c disc_counts_filtered/chrom/sample_name.txt.gz) <(bgzip -d -c tloc_counts_filtered/chrom/sample_name.txt.gz) | bgzip -c > sample_counts/chrom/sample_name.txt.gz
```

