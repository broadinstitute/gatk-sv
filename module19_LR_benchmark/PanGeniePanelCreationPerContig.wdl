version 1.0

import "Structs.wdl"


workflow PanGeniePanelCreationPerContig {
    input {
        File phased_bcf
        File phased_bcf_idx
        File reference_fasta
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing = 0.2
        String output_prefix

        String docker
        File? monitoring_script
    }

    call PanGeniePanelCreation {
        input:
            phased_bcf = phased_bcf,
            phased_bcf_idx = phased_bcf_idx,
            reference_fasta = reference_fasta,
            prepare_vcf_script = prepare_vcf_script,
            add_ids_script = add_ids_script,
            merge_vcfs_script = merge_vcfs_script,
            frac_missing = frac_missing,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script
    }

    output {
        File panel_vcf_gz = PanGeniePanelCreation.panel_vcf_gz
        File panel_vcf_gz_tbi = PanGeniePanelCreation.panel_vcf_gz_tbi
        File panel_id_split_vcf_gz = PanGeniePanelCreation.panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi = PanGeniePanelCreation.panel_id_split_vcf_gz_tbi
    }
}


# TODO consider piping more steps
task PanGeniePanelCreation {
    input {
        File phased_bcf
        File phased_bcf_idx
        File reference_fasta
        String output_prefix

        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing

        String docker
        File? monitoring_script

        RuntimeAttr? runtime_attr_override
    }

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # validate variants against reference
        bcftools norm --check-ref e --fasta-ref ~{reference_fasta} ~{phased_bcf} &> validate-vcf.log

        # run PanGenie prepare-vcf script
        bcftools view ~{phased_bcf} | \
            python3 ~{prepare_vcf_script} \
                --missing ~{frac_missing} \
            2> prepare-vcf.log \
            1> prepare.vcf

        bcftools stats prepare.vcf > ~{output_prefix}.prepare.stats.txt

        # run PanGenie add-ids script
        cat prepare.vcf | \
            python3 ~{add_ids_script} \
            2> add-ids.log \
            1> prepare.id.vcf

        # split to biallelic
        bcftools norm -m- prepare.id.vcf \
            --threads $(nproc) \
            2> split.log \
            1> prepare.id.split.vcf

        # run PanGenie merge script
        pip install pyfaidx
        python3 ~{merge_vcfs_script} merge \
            -vcf prepare.id.split.vcf \
            -r ~{reference_fasta} \
            -ploidy 2  \
            2> merge-haplotypes.log \
            1> prepare.id.split.mergehap.vcf

        # PanGenie script emits header with missing contig lines, so we must bgzip and index
        bcftools view prepare.id.split.mergehap.vcf \
            -Oz -o ~{output_prefix}.prepare.id.split.mergehap.vcf.gz
        bcftools index -t ~{output_prefix}.prepare.id.split.mergehap.vcf.gz

        bcftools stats ~{output_prefix}.prepare.id.split.mergehap.vcf.gz > ~{output_prefix}.prepare.id.split.mergehap.stats.txt

        bcftools sort prepare.id.split.vcf -Oz -o ~{output_prefix}.prepare.id.split.vcf.gz
        bcftools index -t ~{output_prefix}.prepare.id.split.vcf.gz
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(phased_bcf, "GiB")*2),
        disk_gb: 15 + ceil(size(phased_bcf, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        }
    output {
        File monitoring_log = "monitoring.log"
        Array[File] logs = glob("*.log")
        File prepare_stats = "~{output_prefix}.prepare.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf_gz = "~{output_prefix}.prepare.id.split.mergehap.vcf.gz"
        File panel_vcf_gz_tbi = "~{output_prefix}.prepare.id.split.mergehap.vcf.gz.tbi"
        File panel_id_split_vcf_gz = "~{output_prefix}.prepare.id.split.vcf.gz"
        File panel_id_split_vcf_gz_tbi = "~{output_prefix}.prepare.id.split.vcf.gz.tbi"
    }
}