# SrTest
This repository contains the workflow that evaluates split read support for all SV calls on a per-batch basis.

## Required matrices
Sr matrices should be prepared for this process. The matrix describes soft clipped alignments in all individuals, and can be collected from the aligned sequences by following these steps:

1. run `svtk collect-pesr` to collect split read information:
```
svtk collect-pesr sample.bam split_count/sample.txt pe_count/sample.txt
```

2. add sample name as an extra column to each pe_count output:
```
python script/add_sample_name.py split_count/sample.txt
```

Now the split_count/sample.txt looks like this:
```
1	9997	left	1	sample
1	9999	left	3	sample
1	10000	left	1	sample
1	10001	left	5	sample
1	10002	left	9	sample
1	10003	left	10	sample
```

3. concatenate the split_count output of all samples, sort, bgzip and tabix:
```
cat split_count/*.txt > matircs.sr.txt
sort -k1,1 -k2,2n  matircs.sr.txt >  matircs.sr.sorted.txt
bgzip -c matircs.sr.sorted.txt > matircs.sr.sorted.txt.gz
tabix -b 2 -e 2 matircs.sr.sorted.txt.gz
```

## Input files
Input files for PeTest are produced through the 01_algorithm_integration step, and are usually kept under `../../01_algorithm_integration/vcfcluster/`. Names of input files are in this format: `{batch}.{source}.{chrom}.vcf.gz`


## Manual Process
Autosomes and allosomes should be processed separately, with two whitelists contaning samples names of all males (whitelists/{batch}.males.list) and females(whitelists/{batch}.females.list) prepared. The whitelists have one sample name in each line.
For autosomes:
```
svtk pe-test ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz matircs.sr.sorted.txt.gz srtest/{batch}.{source}.{chrom}.stats
```
For allosomes:
```
svtk sr-test --samples whitelists/{batch}.females.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.sr.sorted.txt.gz srtest_allosomes/{batch}.{source}.{chrom}.females.stats
svtk sr-test --samples whitelists/{batch}.males.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.sr.sorted.txt.gz srtest_allosomes/{batch}.{source}.{chrom}.males.stats
python script/sr_merge_allosomes.py batch source chrom X
python script/sr_merge_allosomes.py batch source chrom Y
```

## Efficient manual process
It is recommended that big vcf input be split randomly into smaller files (e.g. 100 SV records per vcf), run through sr-test and then merge the splits together. To split vcf files:
```
python script/split_vcf.random.py ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz split_vcf/{batch}.{source}.{chrom}. -s number_of_svs_per_split
```
Sort and index each split vcf (i.e. `split_out`),
```
vcf-sort split_out > split_out.vcf
bgzip split_out.vcf
tabix split_out.vcf.gz
```
For each split vcf, apply `svtk sr-test` and them merge them:
```
svtk sr-test split_out.vcf.gz matircs.split.sorted.txt.gz split_srtest/split_out.stats
cat {input} | sed -r -e '/^chr\\s/d' | sort -k1,1V -k2,2n | cat <(head -n1 {input[0]}) - > {output}
```


Here's full instruction of `script/split_vcf.random.py`:
```
usage: split_vcf.random.py [-h] [-s SIZE] input output

positional arguments:
  input                 namd of input vcf.gz to be splited
  output                prefix of output

optional arguments:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  size of outputs
```

## Output filesthi sformat:
Result from this step are kept under the `srtest/` folder with names in the format: `{batch}.{source}.{chrom}.stats`

