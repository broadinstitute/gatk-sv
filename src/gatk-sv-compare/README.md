**THIS IS A ROUGH DRAFT**

A new package for directly comparing gatk-sv VCFs. It should be written as a python package called "gatk-sv-compare". The goal is to be able to compare any two VCFs from gatk-sv and for an analyst to be able to easily understand how they are similar and how they are different. 

Look at the MainVcfQc.wdl for examples of methods in QC. I want the following to compose the analysis. 

- AF=allele frequency, can be bucketed into singleton, AC >1 and AF < 0.01, and AF > 0.01 when bucketing AF.
- Genomic context will be assigned using bed tracks with the gatk tool at /Users/markw/IdeaProjects/gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/sv/SVRegionOverlap.java - these tracks are to be provided as --seg-dup-track, --simple-repeat-track, --repeatmasker-track bed files. 
- SV size buckets should only apply to DEL and DUP types with ranges 50-500, 500-5000, and >5000 bp. 
- Mobile element insertions (MEIs) with svtype INS but alt <SV:ME:Alu> for example should grouped into a single "INS:MEI" svtype.

The first subcommand "preprocess" takes in the two vcfs and preprocesses them with gatk using SVConcordance: two runs, each swapping the --eval and --truth vcf inputs. It then runs each output concordance annotated vcf through SVRegionOverlap with the provided tracks.

Then run each of the vcf outputs from SVRegionOverlap through each of the following analyses and put outputs in corresponding subdirectories in a specified output directory. All intermediate files should be put into the output directory as well.

- Site-level comparisons: 
  - binned_counts: Compressed tables of overall counts by svtype x svsize x genomic_context x AF bucket. i.e. a table for small DEL in simple repeats with AF<0.01, etc.
  - overall_counts: Visualization of svtype count, svsize distributions by svtype, AF distributions by svtype, and genomic context by svtype in each call set
  - genotype_dist: Visualization of hardy-weinberg distributions as ternary plots. Statistics of number of sites in HWE nominally and after bonferroni correction. 
  - counts_per_genome: Visualization of SV counts per genome
  - site_overlap_plots: Bar charts depicting counts of overlapping sites and mutually exclusive sites (according to SVConcordance STATUS info field) by overall svtype, svize, genomic context, and AF bucket. 
  - site_overlap_tables: Tables of overlap metrics for all svtype x svsize x genomic_context x AF bucket combinations.
  - allele_freq: allele frequency correlations - plotted overall and for each svtype, svize, genomic context, and AF bucket. Tables for svtype x svsize x genomic_context x AF bucket combinations.
  - genotype_quality : histograms of genotype quality for each svtype, svize, genomic context, and AF bucket. Tables for svtype x svsize x genomic_context x AF bucket combinations.

- If the VCFs contain any overlapping samples, also produce a site-level analysis:
  - genotype_exact_match: Distributions of genotype exact matching by svtype/svsize/genomic_context/AF_bucket, each separately, visualized as histograms. Tables of genotype exact matching by svtype x svsize x genomic_context x AF bucket combinations.
    - Repeat this for not just exact genotype match, but hom ref <-> het changes, hom alt <-> hom ref changes, and het <-> hom ref changes each as well.
  - genotype_concordance_metrics: Distributions of SVConcordance INFO fields including VAR_SENSITIVITY, VAR_PPV, HET_SENSITIVITY, GENOTYPE_CONCORDANCE, etc. as histograms for each svtype, svize, genomic context, and AF bucket


use the following python libraries:
- pysam
- matplotlib
- numpy
- pandas

data input info:
- vcfs may be very large, containing 100k samples and/or ~10M variants. 

optimizations:
The tools must run as fast as possible. Do not duplicate work and do not duplicate logic (DRY). DO NOT try to load the entire vcfs into memory. Break each task down into chunks - it can often be useful to chunk by chromosome. Do not load huge pandas dataframes. When constructing pandas dataframes, do so efficiently and not with repeated insertions. Parallelize where efficient, as this will often be run on multicore machines.
