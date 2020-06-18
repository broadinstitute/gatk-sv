# PeTest
This repository contains the workflow that evaluates pair-end support for all SV calls on a per-batch basis.

## Required matrics
Pe matrics should be prepared for this process. The matircs describes discordant read pairs in all individuals, and can be collectd from the aligned sequenes by following these steps:

1. run `svtk collect-pesr` to collect discordant pair end information:	
```
svtk collect-pesr sample.bam split_count/sample.txt pe_count/sample.txt
```

2. add sample name as an extra column to each pe_count output:
```
python script/add_sample_name.py pe_count/sample.txt
```

Now the pe_count/sample.txt looks like this:
```
1	9995	-	1	249240383	-	sample
1	9997	-	3	197900294	-	sample
1	9999	-	1	10179	+	sample
1	9999	-	5	11254	+	sample
```
3. concatenate the pe_count output of all samples, sort, bgzip and tabix:
```
cat pe_count/*.txt > matircs.pe.txt
sort -k1,1 -k2,2n  matircs.pe.txt >  matircs.pe.sorted.txt
bgzip -c matircs.pe.sorted.txt > matircs.pe.sorted.txt.gz
tabix -b 2 -e 2 matircs.pe.sorted.txt.gz
```

## Input files
Input files for PeTest are produced through the 01_algorithm_integration step, and are usually kept under `../../01_algorithm_integration/vcfcluster/`. Names of input files are in this format: `{batch}.{source}.{chrom}.vcf.gz`


## Manual process
Autosomes and allosomes should be processed separately, with two whitelists contaning samples names of all males (whitelists/{batch}.males.list) and females(whitelists/{batch}.females.list) prepared. The whitelists have one sample name in each line. 

For autosomes:
```
svtk pe-test ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz matircs.pe.sorted.txt.gz petest/{batch}.{source}.{chrom}.stats
```
For allosomes:
```
svtk pe-test --samples whitelists/{batch}.females.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.females.stats
svtk pe-test --samples whitelists/{batch}.males.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.males.stats
python script/pe_merge_allosomes.py batch source chrom X
python script/pe_merge_allosomes.py batch source chrom Y
```

## Efficient manual process
It is recommended that big vcf input be split randomly into smaller files (e.g. 100 SV records per vcf), run through pe-test and then merge the splits together. To split vcf files:
```
python script/split_vcf.random.py ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz split_vcf/{batch}.{source}.{chrom}. -s number_of_svs_per_split
```
Sort and index each split vcf (i.e. `split_out`),
```
vcf-sort split_out > split_out.vcf
bgzip split_out.vcf
tabix split_out.vcf.gz
```
For each split vcf, apply `svtk pe-test` and them merge them:
```
svtk pe-test split_out.vcf.gz matircs.pe.sorted.txt.gz split_petest/split_out.stats
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

## Output file format: 
Result from this step are kept under the `petest/` folder with names in the format: `{batch}.{source}.{chrom}.stats` 



