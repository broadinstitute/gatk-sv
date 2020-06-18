version 1.0

import "Structs.wdl"

workflow make_cohort_VCFs {
  input{
    Array[File] pesr_vcfs
    Array[File] depth_vcfs
    String sv_pipeline_docker
  }
  call MergePESRVcfs {
    input:
      vcfs_list=write_lines(pesr_vcfs),
      prefix="pesr",
      sv_pipeline_docker=sv_pipeline_docker
  }
  call MergeDepthVcfs{
    input:
      vcfs=depth_vcfs,
      prefix="depth",
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File pesrlookup = MergePESRVcfs.lookup
    File depthlookup = MergeDepthVcfs.lookup
  }
}


task MergePESRVcfs {
  input{
    File vcfs_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 64,
    disk_gb: 200,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.sh ~{vcfs_list} ~{prefix}
    while read vcf; do
        svtk vcf2bed --no-header $vcf $vcf.bed
        awk '{if($3-$2>5000) print $0}' $vcf.bed > test.bed;mv test.bed $vcf.bed
    done < vcfs.list
    cat *.bed |sort -k1,1V -k2,2n -k3,3n |bgzip -c > cohort.sort.bed.gz
    #clustering
    svtk vcf2bed ~{prefix}.vcf.gz ~{prefix}.vcf.gz.bed
    awk '{if($3-$2>5000) print $0}' ~{prefix}.vcf.gz.bed >test.bed ; mv test.bed ~{prefix}.vcf.gz.bed
    fgrep DEL ~{prefix}.vcf.gz.bed> del.bed
    fgrep DUP ~{prefix}.vcf.gz.bed> dup.bed
    svtk bedcluster del.bed |cut -f1-7 |awk '{print $0","}' > del.cluster.bed
    svtk bedcluster dup.bed |cut -f1-7 |awk '{print $0","}' > dup.cluster.bed
    cat del.cluster.bed dup.cluster.bed |sort -k1,1V -k2,2n -k3,3n |fgrep -v "#"> cluster.bed
    #harrison's 
    zcat cohort.sort.bed.gz | awk '{a[$1"@"$2"@"$3]=a[$1"@"$2"@"$3]?a[$1"@"$2"@"$3]":"$4:$4;b[$1"@"$2"@"$3]=b[$1"@"$2"@"$3]?b[$1"@"$2"@"$3]","$6:$6;}END{for (i in a)print i "\t" a[i] "\t" b[i];}'|tr '@' '\t'|bgzip>pesr.combined.gz
    zcat pesr.combined.gz |awk -F "," ' { print $0"\t"NF } ' |awk -v OFS="\t" '{print $1,$2,$3,$4":",$6}' |bgzip -c > pesr.lookup.gz
    >>>
  output {
    File cohort_sort = "cohort.sort.bed.gz"
    File cluster="cluster.bed"
    File lookup="pesr.lookup.gz"
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

task MergeDepthVcfs {
  input {
    Array[File] vcfs
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py ~{write_lines(vcfs)} ~{prefix}.vcf
    vcf-sort -c ~{prefix}.vcf | bgzip -c > ~{prefix}.vcf.gz
    while read vcf; do
        local_vcf=$(basename $vcf)
        svtk vcf2bed --no-header $vcf $local_vcf.bed   # for each depth vcf make bed, duplicated
    done < ~{write_lines(vcfs)}
    cat *.bed |sort -k1,1V -k2,2n -k3,3n> cohort.sort.bed # concat raw depth vcf, duplicated
    svtk vcf2bed ~{prefix}.vcf.gz ~{prefix}.vcf.gz.bed   # vcf2bed merge_vcfs, non_duplicated
    fgrep DEL ~{prefix}.vcf.gz.bed> del.bed  # del non duplicated
    fgrep DUP ~{prefix}.vcf.gz.bed> dup.bed  # dup non duplicated
    svtk bedcluster del.bed |cut -f1-7 |awk '{print $0","}' > del.cluster.bed #cluster non_duplicated del
    svtk bedcluster dup.bed |cut -f1-7 |awk '{print $0","}' > dup.cluster.bed #cluster non_duplicated dup
    cat del.cluster.bed dup.cluster.bed |sort -k1,1V -k2,2n -k3,3n |fgrep -v "#"> cluster.bed #combine clusterd non-duplicated
    python3 <<CODE
    varID={}
    with open("cohort.sort.bed",'r') as f: # From the depth cohort bed, a dictionary of all duplicate variants and their samples
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
    with open("cohort.combined.bed",'w') as f: # For each unique variant a line with variants and samples
        for variant in varID.keys():
            CHROM=variant.split(":")[0]
            START=variant.split(':')[1].split("-")[0]
            END=variant.split(':')[1].split("-")[1].split('_')[0]
            varcol=":".join(varID[variant]["varids"])+':'
            samplecol=",".join(varID[variant]['sample'])+','
            f.write(CHROM+"\t"+START+"\t"+END+"\t"+varcol+"\t"+samplecol+'\n')
    samp={}
    var={}
    with open("cohort.combined.bed","r") as f: # for EACH variant ID, a list of duplicate variants and samples
        for line in f:
    #         print(line)
            dat=line.split('\t')        
            for variant in dat[3].split(":")[0:-1]:
                samp[variant]=dat[4].split(',')[0:-1]
                var[variant]=dat[3].split(":")[0:-1]
    with open("master_cluster_dups.bed",'w') as g: # Using cluster.bed, for each 
        with open("cluster.bed","r") as f:
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
    with open("cluster.combined.bed",'w') as g:
        with open("cluster.bed","r") as f:
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
    bgzip master_cluster_dups.bed
  >>>
  output {
    File lookup = "master_cluster_dups.bed.gz"
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
