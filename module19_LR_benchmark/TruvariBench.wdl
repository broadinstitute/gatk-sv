version 1.0
import "Structs.wdl"

workflow CallTruvariBench {
    input {
        Array[String] chromosomes
        File truth_vcf
        File comp_vcf
        String prefix
        File ref_dict

        String? truvari_params
    }

    scatter (chr in chromosomes) {
        call TruvariBench {
            input:
                chromosome = chr,
                truth_vcf = truth_vcf,
                comp_vcf = comp_vcf,
                params = truvari_params
        }
    }

    call MergePerChrCalls as merge_tp_base {
     input: 
       vcfs=TruvariBench.tp_base_vcf, 
       prefix="~{prefix}.tp_base", 
       ref_dict=ref_dict 
     }

    call MergePerChrCalls as merge_tp_comp {
     input: 
      vcfs=TruvariBench.tp_comp_vcf, 
      prefix="~{prefix}.tp_comp", 
      ref_dict=ref_dict 
      }

    call MergePerChrCalls as merge_fp { 
    input: 
    vcfs=TruvariBench.fp_vcf, 
    prefix="~{prefix}.fp", 
    ref_dict=ref_dict }

    call MergePerChrCalls as merge_fn { 
    input: 
    vcfs=TruvariBench.fn_vcf, 
    prefix="~{prefix}.fn", 
    ref_dict=ref_dict 
    }
    
    output {
        File tp_base_vcf = merge_tp_base.merged_vcf
        File tp_comp_vcf = merge_tp_comp.merged_vcf
        File fp_vcf = merge_fp.merged_vcf
        File fn_vcf = merge_fn.merged_vcf
    }
}

task TruvariBench {
    input {
        String chromosome
        File truth_vcf
        File comp_vcf
        String params=""

        RuntimeAttr? runtime_attr_override
    }
   
    Int disk_size = ceil(size(truth_vcf, "GB") + size(comp_vcf, "GB")) * 3 + 50

    command <<<
        set -euxo pipefail
        
        tabix ~{truth_vcf}
        tabix ~{comp_vcf}

        bcftools view -Oz ~{truth_vcf} ~{chromosome} > truth_~{chromosome}.vcf.gz
        tabix truth_~{chromosome}.vcf.gz

        bcftools view -Oz ~{comp_vcf} ~{chromosome} > comp_~{chromosome}.vcf.gz
        tabix comp_~{chromosome}.vcf.gz

        truvari bench -b truth_~{chromosome}.vcf.gz -c comp_~{chromosome}.vcf.gz -o outdir ~{params}
    >>>
    
    output {
        File tp_base_vcf = "outdir/tp-base.vcf.gz"
        File tp_comp_vcf = "outdir/tp-comp.vcf.gz"
        File fp_vcf = "outdir/fp.vcf.gz"
        File fn_vcf = "outdir/fn.vcf.gz"
    }

     #########################
     RuntimeAttr default_attr = object {
         cpu_cores:          4,
         mem_gb:             8,
         disk_gb:            disk_size,
         boot_disk_gb:       10,
         preemptible_tries:  2,
         max_retries:        0,
         docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
     }
     RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
     runtime {
         cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
         memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
         disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
         bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
         preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
         maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
         docker:                 select_first([runtime_attr.docker,            default_attr.docker])
     }
}



task MergePerChrCalls {
  input {
    Array[File] vcfs  # Input list of bgzipped VCFs (vcf.gz)
    File ref_dict
    String prefix

    RuntimeAttr? runtime_attr_override
  }

   Int disk_size = 2*ceil(size(vcfs, "GB")) + 10

  command <<<
    # Index each VCF using tabix
    for vcf in ~{sep=' ' vcfs}; do
      tabix -p vcf ${vcf}
    done

    # Concatenate using bcftools
    bcftools concat -a -Oz -o ~{prefix}.vcf.gz ~{sep=' ' vcfs}
    tabix -p vcf ~{prefix}.vcf.gz
  >>>

  output {
    File merged_vcf = "~{prefix}.vcf.gz"
    File merged_vcf_index = "~{prefix}.vcf.gz.tbi"
  }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
