#workdir: /data/talkowski/xuefang/data/gnomad_V3/module07/cleanvcf/per_chr_info
#mkdir US_RM_SD_SR
awk '{print $1,$2,$2,$4,$5}'  gnomad-sv-v3.chr18.final_cleanup.info3  | sed -e 's/ /\t/g' > US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp
awk '{print $1,$3,$3,$4,$5}'  gnomad-sv-v3.chr18.final_cleanup.info3  | sed -e 's/ /\t/g' > US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp
awk '{if ($8=="DEL" || $8=="DUP" || $8=="CNV" ) print}' gnomad-sv-v3.chr18.final_cleanup.info3 | awk '{if ($3-$2>5000) print}' | cut -f1-5 >  US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv

bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SimpRep.sorted.merged.bed | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp.vs.SR
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SegDup.sorted.merged.bed  | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp.vs.SD
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.RM.sorted.merged.bed  | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.le_bp.vs.RM

bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SimpRep.sorted.merged.bed | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp.vs.SR
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SegDup.sorted.merged.bed  | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp.vs.SD
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.RM.sorted.merged.bed  | awk '{if ($9>0) print}'> US_RM_SD_SR/gnomad-sv-v3.chr18.ri_bp.vs.RM

bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SimpRep.sorted.merged.bed  > US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv.vs.SR
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.SegDup.sorted.merged.bed  > US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv.vs.SD
bedtools coverage -a US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv -b /data/talkowski/xuefang/data/reference/RM_SD/hg38/hg38.RM.sorted.merged.bed  > US_RM_SD_SR/gnomad-sv-v3.chr18.lg_cnv.vs.RM


