version 1.0

import "Genotype_2.wdl" as g2
import "CombineReassess.wdl" as creassess
import "Utils.wdl" as util

workflow RegenotypeCNVs {
  input {
    String sv_base_mini_docker
    String sv_pipeline_docker
    Array[File] depth_vcfs
    File cohort_depth_vcf
    Array[File] batch_depth_vcfs
    Array[File] coveragefiles
    Array[File] coveragefile_idxs
    Array[File] medianfiles
    Array[File] RD_depth_sepcutoffs
    Int n_per_split
    Int n_RdTest_bins
    Array[String] batches
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    File contig_list
    Array[File] regeno_coverage_medians # one file per batch
    Float regeno_max_allele_freq = 0.01 # Rare variant filter for regenotyping candidates: must be < AF threshold (this parameter) or <= AC threshold (below)
    Int regeno_allele_count_threshold = 3 # Rare variant filter for regenotyping candidates: must be < AF threshold (above) or <= AC threshold (this parameter)
    Int min_var_per_sample_outlier_threshold = 3 # Threshold below which regeno SV count per sample should not be considered an outlier (need when counts are sparse)
    Float regeno_sample_overlap = 0.7 # Minimum sample overlap required between raw and regenotyped calls

    RuntimeAttr? runtime_attr_cluster_merged_depth_beds
    RuntimeAttr? runtime_attr_regeno_raw_combined_depth
    RuntimeAttr? runtime_attr_regeno_merged_depth
    RuntimeAttr? runtime_attr_sample_lookup
    RuntimeAttr? runtime_attr_concat_samplecountlookup
    RuntimeAttr? runtime_attr_concat_sampleidlookup

    RuntimeAttr? runtime_attr_merge_list
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_get_count_cohort_samplelist
    RuntimeAttr? runtime_attr_get_regeno
    RuntimeAttr? runtime_attr_get_median_subset
    RuntimeAttr? runtime_attr_median_intersect
    RuntimeAttr? runtime_attr_concat_regenotyped_vcfs

    # Genotype_2
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_add_batch_samples
    RuntimeAttr? runtime_attr_get_regeno_g2
    RuntimeAttr? runtime_attr_split_beds
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_rd_test_gt_regeno
    RuntimeAttr? runtime_attr_integrate_depth_gq
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_concat_regenotyped_vcfs_g2

    #CombineReassess
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_merge_list_creassess 
  }

  call ClusterMergedDepthBeds {
    input:
      cohort_depth_vcf = cohort_depth_vcf,
      cohort = cohort,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_cluster_merged_depth_beds
  }
  
  call MakeRawCombinedBed {
    input:
      vcfs = batch_depth_vcfs,
      cohort = cohort,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_regeno_raw_combined_depth
  }

  call MakeMergedDepthBeds {
    input:
      regeno_raw_combined_depth = MakeRawCombinedBed.regeno_raw_combined_depth,
      cohort = cohort,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_regeno_merged_depth
  }

  Array[Array[String]] contigs = read_tsv(contig_list)
  scatter (contig in contigs) {
    call MakeSampleLookupBeds {
      input:
        regeno_merged_depth = MakeMergedDepthBeds.regeno_merged_depth,
        regeno_merged_depth_clustered = ClusterMergedDepthBeds.regeno_merged_depth_clustered,
        cohort = cohort,
        contig = contig[0],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_sample_lookup
    }
  }

  call ConcatBed as ConcatSampleCountLookupBed {
    input: 
      bed_shards = MakeSampleLookupBeds.regeno_sample_counts_lookup,
      filename = cohort + ".regeno.sample_counts_lookup.bed",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_samplecountlookup
  }

  call ConcatBed as ConcatSampleIdLookupBed {
    input: 
      bed_shards = MakeSampleLookupBeds.regeno_sample_ids_lookup,
      filename = cohort + ".regeno.sample_ids_lookup.bed",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_sampleidlookup
  }

  scatter(i in range(length(batches))) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = depth_vcfs[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_ids_from_vcf
    }
    File depth_vcf_indexes_ = depth_vcfs[i] + ".tbi"
  }

  call GetAndCountCohortSampleList {
    input:
      batch_sample_lists = GetSampleIdsFromVcf.out_file,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_get_count_cohort_samplelist
  }
  scatter(i in range(length(batches))) {
    call GetRegenotype {
      input:
        depth_genotyped_vcf = depth_vcfs[i],
        Batch = batches[i],
        regeno_sample_counts_lookup = ConcatSampleCountLookupBed.concat_bed,
        regeno_raw_combined_depth = MakeRawCombinedBed.regeno_raw_combined_depth,
        n_samples_cohort = GetAndCountCohortSampleList.n_samples_cohort,
        regeno_max_allele_freq = regeno_max_allele_freq, 
        regeno_allele_count_threshold = regeno_allele_count_threshold, 
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_get_regeno
    }
    call GetMedianSubset {
      input: 
        medians = regeno_coverage_medians[i],
        batch = batches[i],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_get_median_subset
    }
    call MedianIntersect {
      input: 
        median = GetMedianSubset.regeno_median,
        regeno_list = GetRegenotype.regeno_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_median_intersect
    }
  }
  call MergeList{
    input:
      regeno_beds = MedianIntersect.regeno_bed,
      prefix="master_regeno",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_merge_list
  }

  if (MergeList.num_regeno > 0) {
    scatter (i in range(length(batches))) {
      call g2.Regenotype as Genotype_2 {
        input:
          depth_vcf=depth_vcfs[i],
          regeno_bed= MergeList.master_regeno,
          cohort_depth_vcf=cohort_depth_vcf,
          batch_depth_vcf=batch_depth_vcfs[i],
          coveragefile=coveragefiles[i],
          coveragefile_idx=coveragefile_idxs[i],
          medianfile=medianfiles[i],
          RD_depth_sepcutoff=RD_depth_sepcutoffs[i],
          n_per_split=n_per_split,
          n_RdTest_bins=n_RdTest_bins,
          batch=batches[i],
          samples_list=GetSampleIdsFromVcf.out_file[i],
          sv_pipeline_docker=sv_pipeline_docker,
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_add_batch_samples = runtime_attr_add_batch_samples,
          runtime_attr_get_regeno_g2 = runtime_attr_get_regeno_g2,
          runtime_attr_split_beds = runtime_attr_split_beds,
          runtime_attr_make_subset_vcf = runtime_attr_make_subset_vcf,
          runtime_attr_rd_test_gt_regeno = runtime_attr_rd_test_gt_regeno,
          runtime_attr_integrate_depth_gq = runtime_attr_integrate_depth_gq,
          runtime_attr_add_genotypes = runtime_attr_add_genotypes,
          runtime_attr_concat_regenotyped_vcfs_g2 = runtime_attr_concat_regenotyped_vcfs_g2
      }
    }

    call creassess.CombineReassess as CombineReassess {
      input:
        samplelist = GetAndCountCohortSampleList.cohort_samplelist,
        regeno_file = MergeList.master_regeno,
        regeno_sample_ids_lookup = ConcatSampleIdLookupBed.concat_bed,
        vcfs = Genotype_2.genotyped_vcf,
        min_var_per_sample_outlier_threshold = min_var_per_sample_outlier_threshold,
        regeno_sample_overlap = regeno_sample_overlap,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_merge_list_creassess = runtime_attr_merge_list_creassess,
        runtime_attr_vcf2bed = runtime_attr_vcf2bed
    }
    
    if (CombineReassess.num_regeno_filtered > 0) {
      scatter (i in range(length(Genotype_2.genotyped_vcf))) {
        call ConcatRegenotypedVcfs {
          input:
            depth_vcf=depth_vcfs[i],
            batch = batches[i],
            regeno_vcf = Genotype_2.genotyped_vcf[i],
            regeno_variants = CombineReassess.regeno_variants,
            runtime_attr_override = runtime_attr_concat_regenotyped_vcfs,
            sv_pipeline_docker = sv_pipeline_docker
        }
      }
    }
  }
  output {
    Array[File] regenotyped_depth_vcfs = select_first([ConcatRegenotypedVcfs.genotyped_vcf, depth_vcfs])
    Array[File] regenotyped_depth_vcf_indexes = select_first([ConcatRegenotypedVcfs.genotyped_vcf_idx, depth_vcf_indexes_])
    File number_regenotyped_file = MergeList.num_regeno_file
    File number_regenotyped_filtered_file = select_first([CombineReassess.num_regeno_filtered_file, MergeList.num_regeno_file])
  }
}

task ClusterMergedDepthBeds {
  input {
    File cohort_depth_vcf
    String cohort
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File regeno_merged_depth_clustered = "~{cohort}.regeno.merged_depth_clustered.bed"
  }
  command <<<
    set -euo pipefail
    svtk vcf2bed ~{cohort_depth_vcf} merged_depth.bed   # vcf2bed merge_vcfs, non_duplicated
    # split DELs and DUPs into separate, non-duplicated BED files. SVTYPE is 5th column of BED
    awk -F "\t" -v OFS="\t" '{ if ($5 == "DEL") { print > "del.bed" } else if ($5 == "DUP") { print > "dup.bed" } }' merged_depth.bed 
    svtk bedcluster del.bed | cut -f1-7 | awk '{print $0","}' > del.cluster.bed #cluster non_duplicated del
    svtk bedcluster dup.bed | cut -f1-7 | awk '{print $0","}' > dup.cluster.bed #cluster non_duplicated dup
    cat del.cluster.bed dup.cluster.bed | sort -k1,1V -k2,2n -k3,3n | fgrep -v "#" > ~{cohort}.regeno.merged_depth_clustered.bed #combine clusterd non-duplicated
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task MakeRawCombinedBed {
  input {
    Array[File] vcfs
    String cohort
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File regeno_raw_combined_depth = "~{cohort}.regeno.raw_combined_depth.bed"
  }
  command <<<
    set -euxo pipefail
    while read vcf; do
        local_vcf=$(basename $vcf)
        svtk vcf2bed --no-header $vcf $local_vcf.bed   # for each depth vcf make bed, duplicated
    done < ~{write_lines(vcfs)}
    rm ~{sep=' ' vcfs}
    cat *.bed | sort -k1,1V -k2,2n -k3,3n > ~{cohort}.regeno.raw_combined_depth.bed # concat raw depth vcf, duplicated
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task MakeMergedDepthBeds {
  input {
    File regeno_raw_combined_depth
    String cohort
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File regeno_merged_depth="~{cohort}.regeno.merged_depth.bed"
  }
  command <<<
    set -euxo pipefail
    python3 <<CODE
    varID={}
    with open("~{regeno_raw_combined_depth}",'r') as f: # From the depth cohort bed, a dictionary of all duplicate variants and their samples
        for line in f:
            dat=line.rstrip().split('\t')
            samples=dat[-1].split(",")
            var=dat[3]
            ID=dat[0]+":"+dat[1]+'-'+dat[2]+'_'+dat[4]
            if ID not in varID.keys():
                varID[ID]={"sample":samples,"varids":[var]}
            else:
                varID[ID]['sample']=varID[ID]['sample']+samples
                varID[ID]['varids'].append(var)
    with open("~{cohort}.regeno.merged_depth.bed",'w') as f: # For each unique variant a line with variants and samples
        for variant in varID.keys():
            CHROM=variant.split(":")[0]
            START=variant.split(':')[1].split("-")[0]
            END=variant.split(':')[1].split("-")[1].split('_')[0]
            varcol=":".join(varID[variant]["varids"])+':'
            samplecol=",".join(varID[variant]['sample'])+','
            f.write(CHROM+"\t"+START+"\t"+END+"\t"+varcol+"\t"+samplecol+'\n')
    CODE
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task MakeSampleLookupBeds {
  input {
    File regeno_merged_depth
    File regeno_merged_depth_clustered
    String cohort
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File regeno_sample_counts_lookup = "~{contig}.regeno.sample_counts_lookup.bed"
    File regeno_sample_ids_lookup = "~{contig}.regeno.sample_ids_lookup.bed"
  }
  command <<<
    set -euxo pipefail
    # select rows of BED files pertaining to contig - chrom is 1st column of each BED file
    touch ~{contig}.regeno.merged_depth.bed ~{contig}.regeno.merged_depth_clustered.bed
    awk -F "\t" -v OFS="\t" '{ if ($1 == "~{contig}") { print > "~{contig}.regeno.merged_depth.bed" } }' ~{regeno_merged_depth}
    awk -F "\t" -v OFS="\t" '{ if ($1 == "~{contig}") { print > "~{contig}.regeno.merged_depth_clustered.bed" } }' ~{regeno_merged_depth_clustered}
    python3 <<CODE
    # dictionary of (samples, varIDs) for de-duplicated variant for EACH varID corresponding to that unique variant
    varID_data = {} 
    with open("~{contig}.regeno.merged_depth.bed","r") as f: # for EACH variant ID, a list of duplicate variants and samples
        for line in f:
            dat=line.split('\t')
            varIDs_list = dat[3].split(":")[0:-1]
            samples_list = dat[4].split(',')[0:-1]
            for varID in varIDs_list:
                varID_data[varID] = (samples_list, varIDs_list)
    with open("~{contig}.regeno.sample_counts_lookup.bed",'w') as g: # Using clustered merged (de-dupped) depth calls, for each clustered variant get varIDs and samples of the component calls
        with open("~{contig}.regeno.merged_depth_clustered.bed","r") as f:
            for line in f:
                samples=[]
                variants=[]
                dat=line.rstrip().split("\t")
                for varID in dat[6][0:-1].split(','):
                    samples.extend(varID_data[varID][0]) # samples are first in tuple
                    samples = list(set(samples))
                    variants.extend(varID_data[varID][1]) # variant IDs are 2nd in tuple
                    variants = list(set(variants))
                g.write(dat[0]+'\t'+dat[1]+'\t'+dat[2]+'\t'+dat[3]+'\t'+dat[4]+'\t'+dat[5]+'\t'+":".join(variants)+':\t'+str(len(samples))+'\n')
    with open("~{contig}.regeno.sample_ids_lookup.bed",'w') as g:
        with open("~{contig}.regeno.merged_depth_clustered.bed","r") as f:
            for line in f:
                samples=[]
                variants=[]
                dat=line.rstrip().split("\t")
                for varID in dat[6][0:-1].split(','):
                    samples.extend(varID_data[varID][0]) # samples are first in tuple
                    samples = list(set(samples))
                    variants.extend(varID_data[varID][1]) # variant IDs are 2nd in tuple
                    variants = list(set(variants))
                g.write(dat[0]+'\t'+dat[1]+'\t'+dat[2]+'\t'+dat[3]+'\t'+dat[4]+'\t'+dat[5]+'\t'+":".join(variants)+':\t'+','.join(samples)+',\n')
    CODE
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ConcatBed {
  input {
    Array[File] bed_shards
    String filename
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File concat_bed = "~{filename}"
  }
  command <<<
    set -euxo pipefail
    while read bed_shard; do
      cat $bed_shard >> ~{filename} # all BED files are headless and sorted so can just concatenate in order
    done < ~{write_lines(bed_shards)}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetAndCountCohortSampleList {
  input {
    Array[File] batch_sample_lists
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    # concatenate batch sample lists (one sample ID per line) into cohort sample list
    while read batch_list; do
      cat $batch_list >> cohort.sample.list
    done < ~{write_lines(batch_sample_lists)}
    # count lines = number of samples
    wc -l < cohort.sample.list | tr -d ' ' > cohort_sample_count.txt
  >>>
  output {
    File cohort_samplelist = "cohort.sample.list"
    File cohort_samplecount_file = "cohort_sample_count.txt"
    Int n_samples_cohort = read_int("cohort_sample_count.txt")
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetRegenotype {
  input {
    File depth_genotyped_vcf
    File regeno_sample_counts_lookup
    File regeno_raw_combined_depth
    Int n_samples_cohort
    Float regeno_max_allele_freq # default = 0.01 set in RegenotypeCNVs.wdl
    Int regeno_allele_count_threshold # default = 3 set in RegenotypeCNVs.wdl
    String Batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail

    # Inputs
    n_samp=~{n_samples_cohort}
    max_af=~{regeno_max_allele_freq}
    min_count=~{regeno_allele_count_threshold}

    # svtk vcf2bed step
    svtk vcf2bed ~{depth_genotyped_vcf} ~{Batch}.bed

    # restrict to variants originating in this batch
    fgrep "~{Batch}_" ~{regeno_raw_combined_depth} > ~{Batch}.origin.raw_combined_depth.bed

    # Python script does the rest
    python3 <<CODE
import sys

batch = "~{Batch}"
n_samples_cohort = int("~{n_samples_cohort}")
max_af = float("~{regeno_max_allele_freq}")
min_count = int("~{regeno_allele_count_threshold}")

bed_file = f"{batch}.bed"
origin_file = f"{batch}.origin.raw_combined_depth.bed"
sample_counts_file = "~{regeno_sample_counts_lookup}"

# ------------------------------------------------------------------------------
# Build maps: identifier -> line fields
def build_variant_map(filepath):
    variant_map = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            identifier = "_".join(fields[0:3] + [fields[4]])
            variant_map[identifier] = fields
    return variant_map

bed_map = build_variant_map(bed_file)
origin_map = build_variant_map(origin_file)

# ------------------------------------------------------------------------------
# Mimic `join`: build lines like `identifier <bed fields> <origin samples>`
joined_lines = []
for var_id in bed_map:
    if var_id in origin_map:
        bed_fields = bed_map[var_id]
        origin_samples = origin_map[var_id][5] if len(origin_map[var_id]) > 5 else ""
        joined_lines.append([var_id] + bed_fields + [origin_samples])

# ------------------------------------------------------------------------------
# Find missing or gained samples
missing_variants = set()
for line in joined_lines:
    line = line + [""] * (8 - len(line))
    identifier, chrom, start, end, var, typ, samples_post, samples_pre = line[:8]

    post_set = set(samples_post.split(",")) if samples_post else set()
    pre_set = set(samples_pre.split(",")) if samples_pre else set()

    missing_in_post = pre_set - post_set
    missing_in_pre = post_set - pre_set

    if int(end) - int(start) > 10000 and typ != "CN0":
        if missing_in_post or missing_in_pre:
            missing_variants.add((chrom, start, end, var, typ))

# ------------------------------------------------------------------------------
# Load sample counts file into lookup
sample_counts = {}
with open(sample_counts_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            var_id = parts[0]
            count = parts[7]
            sample_counts[var_id] = count

# ------------------------------------------------------------------------------
# Get a sample name from the BED
sample = next((fields[5].split(",")[0] for fields in bed_map.values() if len(fields) > 5 and fields[5]), "NA")

# ------------------------------------------------------------------------------
# Apply AF / AC filter and write output
with open(f"{batch}.to_regeno.bed", "w") as out:
    for chrom, start, end, var, typ in sorted(missing_variants):
        var_lookup = f"{var}:"
        num = int(sample_counts.get(var_lookup, 0))
        af = num / n_samples_cohort if n_samples_cohort else 0
        if af < max_af or num <= min_count:
            out.write(f"{chrom}\t{start}\t{end}\t{var}\t{sample}\t{typ}\n")
CODE
>>>

  output {
    File regeno_bed="~{Batch}.to_regeno.bed"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetMedianSubset {
  input {
    String batch
    File medians # this task runs per batch, so there is just one medians file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 <<CODE

    from statistics import median
    import gzip

    # Flexibly open .gz or uncompressed file to read
    # Expect bgzipped BED file but should be able to handle legacy uncompressed files
    # Output file is smaller and deleted as a workflow intermediate so can be uncompressed
    def _open(filename):
      if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
      else:
        return open(filename, 'r')


    with _open("~{medians}") as inp, open("~{batch}_to_regeno.bed", 'w') as outp:
      for line in inp:
        fields = line.strip().split('\t')
        # first 4 fields are variant info (chr, start, end, varID)
        var_info = fields[0:4]
        if any([x == '.' for x in fields[4:]]):
          # skip variants with empty coverage data
          continue
        # else:
        # last len-4 fields are sample coverage medians for each sample in the batch
        med = median([float(x) for x in fields[4:]])
        if (0.85 < med < 0.97) or (1.03 < med < 1.65):
          outp.write('\t'.join(var_info) + "\n")
    CODE
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
  output {
    File regeno_median = "~{batch}_to_regeno.bed"
  }
}

task MedianIntersect {
  input {
    File median
    File regeno_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 2
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    cut -f 4 ~{median} > medmissingID.txt
    { fgrep -w -f medmissingID.txt ~{regeno_list} || true; }|sort -k1,1V -k2,2n -k3,3n -u > median_missing_sorted.bed
  >>>
  output {
    File regeno_bed="median_missing_sorted.bed"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
} 

task MergeList {
  input {
    String prefix
    Array[File] regeno_beds
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail
    # concatenate, remove duplicate variant IDs (ignore other columns), then sort by position
    cat ~{sep=' ' regeno_beds} | sort -u -k4,4 | sort -k1,1V -k2,2n -k3,3n > ~{prefix}.bed
    # count non-empty lines in regeno bed file to determine if empty or not --> proceed with regenotyping or stop here?
    # the OR clause is to ignore return code = 1 because that isn't an error, it just means there were 0 matched lines (but don't ignore real error codes > 1)
    NUM_REGENO=$(grep -c '[^[:space:]]' ~{prefix}.bed || [[ $? == 1 ]] ) 
    echo $NUM_REGENO > regeno_num_lines.txt
  }
  output {
    File master_regeno="~{prefix}.bed"
    Int num_regeno = read_int("regeno_num_lines.txt")
    File num_regeno_file = "regeno_num_lines.txt"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ConcatRegenotypedVcfs {
  input {
    String batch
    File regeno_variants
    File depth_vcf
    File regeno_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    zcat ~{regeno_vcf} |fgrep "#" > head.txt
    zcat ~{regeno_vcf} |fgrep -wf ~{regeno_variants} >body.txt
    cat head.txt body.txt|bgzip -c > regeno.vcf.gz
    zcat ~{depth_vcf} |fgrep -wf ~{regeno_variants} -v |bgzip -c > no_variant.vcf.gz
    vcf-concat regeno.vcf.gz no_variant.vcf.gz \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.depth.regeno_final.vcf.gz
    tabix ~{batch}.depth.regeno_final.vcf.gz
  >>>
  output {
    File genotyped_vcf = "~{batch}.depth.regeno_final.vcf.gz"
    File genotyped_vcf_idx = "~{batch}.depth.regeno_final.vcf.gz.tbi"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
