
##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/00_pesr_preprocessing_MMDLW/15/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"
# Does perlim translocation resolve from raw manta calls
workflow TinyResolve {
  input {
  	String sample
  	File manta_vcf
  	File discfile
    File cytoband
    File mei_bed
    File contigs
    Int min_svsize
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr
  }

  File discfile_idx = discfile + ".tbi"
  File cytoband_idx = cytoband + ".tbi"


    call StandardizeVCFs {
    	input:
	        vcf = manta_vcf,
	        sample_id = sample,
	        caller = "manta",
	        contigs = contigs,
	        min_svsize = min_svsize,
	        sv_pipeline_docker = sv_pipeline_docker,
	        runtime_attr_override = runtime_attr
	  }

    call ResolveManta {
      input:
        vcf = StandardizeVCFs.out,
        sample_id = sample,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband=cytoband,
        cytoband_idx=cytoband_idx,
        discfile=discfile,
        discfile_idx=discfile_idx,
        mei_bed=mei_bed,
        runtime_attr_override=runtime_attr
    }
  }

  output {
    Array[File] tloc_unresolved_vcf = flatten(ResolveManta.unresolved_vcf)
    Array[File] tloc_manta_vcf = flatten(ResolveManta.tloc_vcf)
  }
}

task ResolveManta {
  input {
    File vcf
    String sample_id
    File cytoband_idx
    File discfile
    File discfile_idx
    File cytoband
    File mei_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(discfile,"GiB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10+input_size),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

      tabix -p vcf $vcf
      bash /opt/sv-pipeline/00_preprocessing/scripts/mantatloccheck.sh $vcf $discfile ${sample_id} ~{mei_bed} ~{cytoband}
      bgzip manta.unresolved.vcf
      mv manta.unresolved.vcf.gz unresolved_${sample_id}.manta.complex.vcf.gz
      mv ${sample_id}.manta.complex.vcf.gz tloc_${sample_id}.manta.complex.vcf.gz
    done
  >>>

  output {
  	File tloc_vcf = "tloc_${sample_id}.manta.complex.vcf.gz"
  	File unresolved_vcf = "unresolved_${sample_id}.manta.complex.vcf.gz"
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

task StandardizeVCFs {
  input {
  	File vcf
    String sample_id
    String caller
    File contigs
    Int min_svsize
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    svtk standardize --sample-names ${sample_id} --prefix ~{caller}_${sample_id} --contigs ~{contigs} --min-size ~{min_svsize} $vcf tmp.vcf ~{caller}
    sample_no=`printf %03d $i`
    bcftools sort tmp.vcf -Oz -o std_${sample_no}.~{caller}.${sample_id}.vcf.gz

  >>>

  output {
    File out = glob("std_*.vcf.gz")
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



