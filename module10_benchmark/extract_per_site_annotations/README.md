# this repository provies scripts to extract per site level characters for boost model  from cleanvcf

### Latest updates:
Site-level annotation has been prepared for each chromosome, including allosomes. The data can be found here: gs://talkowski-sv-gnomad-output/zero/CleanVcf/site_level_anno/


**all rest scripts should be run on erisone**

**recommended working directory on erisone: /data/talkowski/xuefang/data/gnomad_V3/module07/cleanvcf/**
## step1: download data:
## cleanvcf:

bash Step0.download_clean_vcf.sh

## sr_background_fail and  sr_bothside_pass SV lists:

bash Step0.download_bg_fail_SVID.sh

bash Step0.download_sr_bothside_SVID.sh

## stitch_fix vcf:

cd stitch_fix

bash Step1.extract_vcf_sites.sh

bash Step2.vcf2bed.sh

## Step2: convert vcf information info bed format 

### note: erisone may not be able to handel large chromosomes such as chr1 and chr2; splitting the vcf before converting is recommended

bash Step2.vcf2bed.sh 

bash  Step3a.bgzip_bed.sh 

bash Step3b.index_bed.sh

# Step3: remove per sample GT from vcf for easier processes:

bash Step4a.extract_vcf_sites.sh

bash Step4b.vcf_sites_2_bed.sh

# Step3a: update the filter columns with bothside pass and high sr background:

Rscript Step5.integrate_sr_bg_and_bothside.R 

bash Step6.split_SVID_filter.sh

# Step3b: annotate SV sites with gnomic context:

mkdir US_RM_SD_SR

bash Step8a.annotate_info3_with_GC.chr10.sh

bash Step8a.annotate_info3_with_GC.chr11.sh

bash Step8a.annotate_info3_with_GC.chr12.sh

bash Step8a.annotate_info3_with_GC.chr13.sh

bash Step8a.annotate_info3_with_GC.chr14.sh

bash Step8a.annotate_info3_with_GC.chr15.sh

bash Step8a.annotate_info3_with_GC.chr16.sh

bash Step8a.annotate_info3_with_GC.chr17.sh

bash Step8a.annotate_info3_with_GC.chr18.sh

bash Step8a.annotate_info3_with_GC.chr19.sh

bash Step8a.annotate_info3_with_GC.chr1.sh

bash Step8a.annotate_info3_with_GC.chr20.sh

bash Step8a.annotate_info3_with_GC.chr21.sh

bash Step8a.annotate_info3_with_GC.chr22.sh

bash Step8a.annotate_info3_with_GC.chr2.sh

bash Step8a.annotate_info3_with_GC.chr3.sh

bash Step8a.annotate_info3_with_GC.chr4.sh

bash Step8a.annotate_info3_with_GC.chr5.sh

bash Step8a.annotate_info3_with_GC.chr6.sh

bash Step8a.annotate_info3_with_GC.chr7.sh

bash Step8a.annotate_info3_with_GC.chr8.sh

bash Step8a.annotate_info3_with_GC.chr9.sh

bash Step8a.annotate_info3_with_GC.chrY.sh

Rscript Step8b.integrate_anno_with_GC.R

bash Step8c.bgzip_info4.sh

#final output:  **gnomad-sv-v3.chrX.final_cleanup.info4.gz**


