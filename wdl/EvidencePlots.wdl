version 1.0

import "Structs.wdl"

workflow IgvEvidence {
    input {
        File ped_file
        File sample_pe_sr
        String family
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_reformat_pe
        RuntimeAttr? runtime_attr_reformat_sr
        RuntimeAttr? runtime_attr_update_pe_sr
    }

    call ProcessEvidence {
        input:
            ped_file = ped_file,
            sample_pe_sr = sample_pe_sr,
            family = family,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_reformat_pe
    }

    output {
        File pe_files = ProcessEvidence.pe_output
        File sr_files = ProcessEvidence.sr_output
        File updated_sample_pe_sr = ProcessEvidence.updated_sample_pe_sr
    }
}

task ProcessEvidence {
    input {
        File ped_file
        File sample_pe_sr
        String family
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([ped_file, sample_pe_sr], "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size * 2),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File pe_output = "pe_evidence.txt"
        File sr_output = "sr_evidence.txt"
        File updated_sample_pe_sr = "updated_sample_pe_sr.txt"
    }

    command <<<
        set -euo pipefail

        # Get family samples
        awk -v fam="~{family}" '$1==fam {print $2}' ~{ped_file} > family_samples.txt

        # Process PE/SR evidence for this family
        echo "sample	pe_file	sr_file" > pe_evidence.txt
        echo "sample	pe_file	sr_file" > sr_evidence.txt
        echo "sample	pe_file	sr_file" > updated_sample_pe_sr.txt

        while read sample; do
            pe_file=$(awk -v s="$sample" '$1==s {print $2}' ~{sample_pe_sr} || echo "")
            sr_file=$(awk -v s="$sample" '$1==s {print $3}' ~{sample_pe_sr} || echo "")
            
            if [[ -n "$pe_file" && -n "$sr_file" ]]; then
                echo "$sample	$pe_file	$sr_file" >> pe_evidence.txt
                echo "$sample	$pe_file	$sr_file" >> sr_evidence.txt
                echo "$sample	$pe_file	$sr_file" >> updated_sample_pe_sr.txt
            fi
        done < family_samples.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
} 