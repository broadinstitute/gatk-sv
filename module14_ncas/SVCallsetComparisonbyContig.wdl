version 1.0

import "Structs.wdl"

workflow LongReadVsShortRead {
    input{
        File long_read_SV_File
        File short_read_SV_File
        Boolean vcf_2_bed_LR
        Boolean vcf_2_bed_SR
        String prefix_LR
        String prefix_SR
        File contig_list 
        File ref_header
        File query_header
        File src_tar
        String sv_base_mini_docker

    }

    if(vcf_2_bed_LR){
        call Vcf2Bed as Vcf2Bed_LR{
            input:
                vcf = long_read_SV_File,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    if(vcf_2_bed_SR){
        call Vcf2Bed as Vcf2Bed_SR{
            input:
                vcf = short_read_SV_File,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    File long_read_SVs = select_first([Vcf2Bed_LR.bed, long_read_SV_File])
    File short_read_SVs = select_first([Vcf2Bed_SR.bed, short_read_SV_File])

    Array[Array[String]] contigs=read_tsv(contig_list)

    scatter(contig in contigs){
        call ExtractQueryRef as extract_query_ref_LR{
            input:
                bed = long_read_SVs,
                contig = contig[0],
                ref_header = ref_header,
                query_header = query_header,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractQueryRef as extract_query_ref_SR{
            input:
                bed = short_read_SVs,
                contig = contig[0],
                ref_header = ref_header,
                query_header = query_header,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVComparison as LR_vs_SR{
            input:
                src_tar = src_tar,
                query = extract_query_ref_LR.query,
                ref = extract_query_ref_SR.ref,
                prefix = "~{prefix_LR}_vs_~{prefix_SR}.~{contig}",
                contig = contig[0],
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVComparison as SR_vs_LR{
            input:
                src_tar = src_tar,
                query = extract_query_ref_SR.query,
                ref = extract_query_ref_LR.ref,
                prefix = "~{prefix_SR}_vs_~{prefix_LR}.~{contig}",
                contig = contig[0],
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    call ConcatComparisons as concat_comparisons_LR_vs_SR{
        input:
            SV_comparison_list = LR_vs_SR.comparison, 
            prefix = "~{prefix_LR}_vs_~{prefix_SR}",
            sv_base_mini_docker = sv_base_mini_docker
        }

    call ConcatComparisons as concat_comparisons_SR_vs_LR{
        input:
            SV_comparison_list = SR_vs_LR.comparison, 
            prefix = "~{prefix_SR}_vs_~{prefix_LR}",
            sv_base_mini_docker = sv_base_mini_docker
        }    


    output{
        File LR_vs_SR_output = concat_comparisons_LR_vs_SR.Concat_file
        File SR_vs_LR_output = concat_comparisons_SR_vs_LR.Concat_file
    }
}


task ConcatComparisons {
    input {
        Array[File] SV_comparison_list
        String prefix
        Boolean? index_output
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    } 
    
    Boolean call_tabix = select_first([index_output, true])
    
    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    Float input_size = size(SV_comparison_list, "GB")
    RuntimeAttr runtime_default = object {
        mem_gb: 2.0,
        disk_gb: ceil(10.0 + input_size * 7.0),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eux

        set -o pipefail

        while read SPLIT; do
          zcat $SPLIT
        done < ~{write_lines(SV_comparison_list)} \
          | (grep -Ev "^#" || printf "") \
          | sort -Vk1,1 -k2,2n -k3,3n \
          | bgzip -c \
          > ~{prefix}.gz

    >>>

  output {
    File Concat_file = "~{prefix}.gz"
  }
}


task Vcf2Bed{
    input{
        File vcf
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed = "~{filebase}.bed.gz"
    }

    String filebase = basename(vcf,".vcf.gz")

    command <<<
        set -Eeuo pipefail

        svtk vcf2bed -i SVTYPE -i SVLEN -i AF --include-filters ~{vcf} ~{filebase}.bed
        bgzip ~{filebase}.bed
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

task ExtractQueryRef{
    input{
        File bed
        File ref_header
        File query_header
        String contig
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File query = "~{contig}.query.gz"
        File ref = "~{contig}.ref.gz"
    }


    command <<<
        set -Eeuo pipefail

        paste <(zcat ~{bed} | cut -f1-3) \
              <(zcat ~{bed}  | cut -f4 | sed -e 's/^/__/' | sed -e 's/$/__/') \
              <(zcat ~{bed} |  cut -f7,8) \
              | awk '{if ($1=="~{contig}") print}' | cat ~{query_header} - | bgzip > ~{contig}.query.gz


        paste <(zcat ~{bed} | cut -f1-3) \
              <(zcat ~{bed} | cut -f4 | sed -e 's/^/__/' | sed -e 's/$/__/') \
              <(zcat ~{bed} | cut -f7,8,9 | sed -e 's/$/\tsample/') \
              | awk '{if ($1=="~{contig}") print}' | cat ~{ref_header} - | bgzip > ~{contig}.ref.gz


   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    }

task SVComparison{
    input{
        File src_tar
        File query
        File ref
        String prefix
        String contig
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File comparison = "~{prefix}.~{contig}.bed.gz"
    }


    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz
        bash src/compare_callsets_V2.sh -O ~{prefix}.~{contig}.bed -p ~{prefix}.~{contig} ~{query} ~{ref} src/
        bgzip ~{prefix}.~{contig}.bed
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
