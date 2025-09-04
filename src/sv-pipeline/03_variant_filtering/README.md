# 03 Variant Filtering
This workflow integrate, filters and genotypes the structural variation(SVs) calls based on the evidence derived in previous modules. The following processes are applied here:
1. Evidences collected from module 02, i.e., rd, pe, sr and baf, are aggregated;
2. Process integrated evidences through Random Forest(RF) to train the optimized parameter for quality control
3. Apply the RF filters on the SVs and remove those failures.

## Manual process
#### Evidence aggragation
a. To aggregate evidence for **pesr callers** (eg. delly, dragen, lumpy, manta, wham), for each `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-r ../02_evidence_assessment/02a_rdtest/rdtest/{batch}.{source}.{chrom}.metrics \
	-p ../02_evidence_assessment/02b_petest/petest/{batch}.{source}.{chrom}.stats \
	-s ../02_evidence_assessment/02c_srtest/srtest/{batch}.{source}.{chrom}.stats \
	-b ../02_evidence_assessment/02d_baftest/baftest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz \
	--fam ../ref/{batch}.fam \
	metrics/{batch}.{source}.{chrom}.metrics
```

b. To aggregate evidence for **read depth callers** (eg. cn.mops, CNVnator, ERDs), for `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-r ../02_evidence_assessment/02a_rdtest/rdtest/{batch}.{source}.{chrom}.metrics \
	-b ../02_evidence_assessment/02d_baftest/baftest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/rdtest_beds/{batch}.{source}.{chrom}.bed \
	--bed \
        --fam ../ref/{batch}.fam \
	metrics/{batch}.{source}.{chrom}.metrics
```
note:  ../ref/{batch}.fam describes family structure of samples included in the study. Here's an example of the format:
`fam_name sample_name father_name mother_name sex case_vs_control`
```
11711	11711.fa	0	0	1	1
11711	11711.mo	0	0	2	1
11711	11711.p1	11711.fa	11711.mo	1	2
11711	11711.s1	11711.fa	11711.mo	1	1
11942	11942.fa	0	0	1	1
11942	11942.mo	0	0	2	1
11942	11942.p1	11942.fa	11942.mo	1	2
11942	11942.s1	11942.fa	11942.mo	1	1
12798	12798.fa	0	0	1	1
12798	12798.mo	0	0	2	1
12798	12798.p1	12798.fa	12798.mo	1	2
12798	12798.s1	12798.fa	12798.mo	2	1
14460	14460.fa	0	0	1	1
14460	14460.mo	0	0	2	1
14460	14460.p1	14460.fa	14460.mo	1	2
14460	14460.s1	14460.fa	14460.mo	2	1
```


c. To aggregate evidence for **mobile element insertion callers** (eg. MELT), for `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-s ../02_evidence_assessment/02c_srtest/srtest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz \
	--batch-list batch_list
	metrics/{batch}.{source}.{chrom}.metrics
```

Modify the aggragated metrics with position of each variants
```
python scripts/add_pos.py -v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz metrics/{batch}.{source}.{chrom}.metrics
python scripts/add_pos.py -v ../01_algorithm_integration/rdtest_beds/{batch}.{source}.{chrom}.bed metrics/{batch}.{source}.{chrom}.metrics
