##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/04_v2_make_cohort_VCFs/3/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"

workflow MergeCohortVcfs {
  input {
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
    RuntimeAttr? runtime_attr_cohort_sort_combined
    RuntimeAttr? runtime_attr_cluster_dups_combined
  }

  call MergeVcfs as MergePESRVcfs {
    input:
      vcfs = pesr_vcfs,
      prefix = cohort + ".all_batches.pesr",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_pesr
  }

  call MergeDepthVcfs {
    input:
      vcfs = depth_vcfs,
      cohort = cohort,
      prefix = cohort + ".all_batches.depth",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_depth
  }

  call MakeCohortSortCombinedBed {
    input:
      vcfs = depth_vcfs,
      cohort = cohort,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_cohort_sort_combined
  }

  call MakeClusterDupsCombinedBed {
    input:
      cohort_combined = MakeCohortSortCombinedBed.cohort_combined,
      cohort_cluster = MergeDepthVcfs.cluster,
      cohort = cohort,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_cluster_dups_combined
  }



  output {
    File cohort_pesr_vcf = MergePESRVcfs.merged_vcf
    File cohort_depth_vcf = MergeDepthVcfs.merged_vcf
    File cohort_combined = MakeCohortSortCombinedBed.cohort_combined
    File cohort_sort = MakeCohortSortCombinedBed.cohort_sort
    File lookup = MakeClusterDupsCombinedBed.lookup
    File cluster_combined = MakeClusterDupsCombinedBed.cluster_combined
  }
}

task MergeVcfs {
  input {
    Array[File] vcfs
    String prefix
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
    File merged_vcf = "~{prefix}.vcf.gz"
  }
  command <<<

    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py ~{write_lines(vcfs)} ~{prefix}.vcf
    rm ~{sep=' ' vcfs}
    vcf-sort -c ~{prefix}.vcf | bgzip -c > ~{prefix}.vcf.gz
  
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

task MergeDepthVcfs {
  input {
    Array[File] vcfs
    String cohort
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged_vcf = "~{prefix}.vcf.gz"
    File cluster="~{cohort}.cluster.bed"
  }
  command <<<

    set -euxo pipefail
    TIMEFORMAT='real: %3lR | user: %3lU | sys: %3lS'
    time /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py ~{write_lines(vcfs)} ~{prefix}.vcf
    time rm ~{sep=' ' vcfs}
    time vcf-sort -c ~{prefix}.vcf | bgzip -c > ~{prefix}.vcf.gz
    time svtk vcf2bed ~{prefix}.vcf.gz ~{prefix}.vcf.gz.bed   # vcf2bed merge_vcfs, non_duplicated
    fgrep DEL ~{prefix}.vcf.gz.bed > del.bed  # del non duplicated
    fgrep DUP ~{prefix}.vcf.gz.bed > dup.bed  # dup non duplicated
    time svtk bedcluster del.bed | cut -f1-7 | awk '{print $0","}' > del.cluster.bed #cluster non_duplicated del
    time svtk bedcluster dup.bed | cut -f1-7 | awk '{print $0","}' > dup.cluster.bed #cluster non_duplicated dup
    cat del.cluster.bed dup.cluster.bed | sort -k1,1V -k2,2n -k3,3n | fgrep -v "#" > ~{cohort}.cluster.bed #combine clusterd non-duplicated
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

task MakeCohortSortCombinedBed {
  input {
    Array[File] vcfs
    String cohort
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File cohort_sort = "~{cohort}.cohort.sort.bed"
    File cohort_combined="~{cohort}.cohort.combined.bed"
  }
  command <<<

    set -euxo pipefail
    TIMEFORMAT='real: %3lR | user: %3lU | sys: %3lS'
    while read vcf; do
        local_vcf=$(basename $vcf)
        svtk vcf2bed --no-header $vcf $local_vcf.bed   # for each depth vcf make bed, duplicated
    done < ~{write_lines(vcfs)}
    rm ~{sep=' ' vcfs}
    cat *.bed | sort -k1,1V -k2,2n -k3,3n > ~{cohort}.cohort.sort.bed # concat raw depth vcf, duplicated
    time python3 <<CODE
    varID={}
    with open("~{cohort}.cohort.sort.bed",'r') as f: # From the depth cohort bed, a dictionary of all duplicate variants and their samples
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
    with open("~{cohort}.cohort.combined.bed",'w') as f: # For each unique variant a line with variants and samples
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

task MakeClusterDupsCombinedBed {
  input {
    File cohort_combined
    File cohort_cluster
    String cohort
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File lookup = "~{cohort}.master_cluster_dups.bed"
    File cluster_combined="~{cohort}.cluster.combined.bed"
  }
  command <<<

    set -euxo pipefail
    TIMEFORMAT='real: %3lR | user: %3lU | sys: %3lS'
    time python3 <<CODE
    samp={}
    var={}
    with open("~{cohort_combined}","r") as f: # for EACH variant ID, a list of duplicate variants and samples
        for line in f:
            dat=line.split('\t')        
            for variant in dat[3].split(":")[0:-1]:
                samp[variant]=dat[4].split(',')[0:-1]
                var[variant]=dat[3].split(":")[0:-1]
    with open("~{cohort}.master_cluster_dups.bed",'w') as g: # Using cluster.bed, for each 
        with open("~{cohort_cluster}","r") as f:
            for line in f:
                samples=[]
                variants=[]
                dat=line.rstrip().split("\t")
                for variant in dat[6][0:-1].split(','):
                    samples=samples+samp[variant]
                    samples=list(set(samples))
                    variants=variants+var[variant]
                    variants=list(set(variants))
                g.write(dat[0]+'\t'+dat[1]+'\t'+dat[2]+'\t'+dat[3]+'\t'+dat[4]+'\t'+dat[5]+'\t'+":".join(variants)+':\t'+str(len(samples))+'\n')
    with open("~{cohort}.cluster.combined.bed",'w') as g:
        with open("~{cohort_cluster}","r") as f:
            for line in f:
                samples=[]
                variants=[]
                dat=line.rstrip().split("\t")
                for variant in dat[6][0:-1].split(','):
                    samples=samples+samp[variant]
                    samples=list(set(samples))
                    variants=variants+var[variant]
                    variants=list(set(variants))
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
