version 1.0
import "Structs.wdl"

workflow MergeVCFLists {
    input {
        Array[File] vcf_list1  # List of bgzipped VCFs
        Array[File] vcf_list2  # List of bgzipped VCFs
        Array[File] vcf_idx_list1  # List of vcf index
        Array[File] vcf_idx_list2  # List of vcf index
        Array[String] contig_list
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_merge1
        RuntimeAttr? runtime_attr_merge2

    }

    # Pairwise merge of VCFs in the two lists
    scatter (i in range(length(vcf_list1))) {
        call ConcatPair {
            input:
                vcf1 = vcf_list1[i],
                vcf2 = vcf_list2[i],
                idx1 = vcf_idx_list1[i],
                idx2 = vcf_idx_list2[i],
                prefix = output_prefix,
                contig = contig_list[i],
                docker_file = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_merge1
        }
    }

    # Merge all pairwise merged VCFs into final output
    call ConcatAll {
        input:
            vcfs_to_merge = ConcatPair.output_vcf,
            idxes_to_merge = ConcatPair.output_idx,
            output_prefix = output_prefix,
            docker_file = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge2

    }

    output {
        File final_vcf = ConcatAll.merged_vcf
        File final_index = ConcatAll.merged_tbi
    }
}

task ConcatPair {
    input {
        File vcf1
        File vcf2
        File idx1
        File idx2
        String prefix
        String contig
        String docker_file
        RuntimeAttr? runtime_attr_override

    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        bcftools concat -a -Oz -o ~{prefix}.~{contig}.vcf.gz ~{vcf1} ~{vcf2}
        #bcftools sort -Oz -o ~{prefix}.~{contig}.vcf.gz ~{prefix}.~{contig}.temp_merged.vcf.gz
        tabix -p vcf ~{prefix}.~{contig}.vcf.gz
        #rm ~{prefix}.~{contig}.temp_merged.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.~{contig}.vcf.gz"
        File output_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task ConcatAll {
    input {
        Array[File] vcfs_to_merge
        Array[File] idxes_to_merge
        String output_prefix
        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<


        for vcf in ~{sep=" " vcfs_to_merge}; do
          echo $vcf >> vcfs.list
        done


        bcftools concat -a -Oz -f vcfs.list -o ~{output_prefix}.vcf.gz
        #bcftools sort -Oz -o ~{output_prefix}.vcf.gz ~{output_prefix}.temp_merged.vcf.gz
        tabix -p vcf ~{output_prefix}.vcf.gz
        #rm ~{output_prefix}.temp_merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_tbi = "~{output_prefix}.vcf.gz.tbi"
    }


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}