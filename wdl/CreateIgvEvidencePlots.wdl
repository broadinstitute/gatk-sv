version 1.0

import "Structs.wdl"
import "EvidencePlots.wdl" as evidence
import "IgvEvidencePlots.wdl" as igv_plots

workflow IgvAllSamples {
    input {
        File ped_file
        File sample_pe_sr
        String buffer
        File? fam_ids
        File varfile
        File reference
        File reference_index
        String prefix
        Boolean is_snv_indel
        Boolean file_localization
        String sv_base_mini_docker
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_run_igv
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_cpx
        RuntimeAttr? runtime_attr_reformat_pe
        RuntimeAttr? runtime_attr_reformat_sr
        RuntimeAttr? runtime_attr_update_pe_sr
    }

    call GetSamples {
        input:
            ped_file = ped_file,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_run_igv
    }

    call GetSamplePeSr {
        input:
            ped_file = ped_file,
            sample_pe_sr = sample_pe_sr,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_run_igv
    }

    scatter (family_name in GetSamples.family_list) {
        call GetFamilySample {
            input:
                ped_file = ped_file,
                family = family_name,
                prefix = prefix,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_override=runtime_attr_run_igv
        }

        call evidence.IgvEvidence as IgvEvidence {
            input:
                ped_file = ped_file,
                sample_pe_sr = GetSamplePeSr.sample_pe_sr_list,
                family = family_name,
                prefix = prefix,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_reformat_pe = runtime_attr_reformat_pe,
                runtime_attr_reformat_sr = runtime_attr_reformat_sr,
                runtime_attr_update_pe_sr = runtime_attr_update_pe_sr
        }

        call igv_plots.Igv as Igv {
            input:
                family = family_name,
                samples = GetFamilySample.family_sample_list,
                varfile = varfile,
                buffer = buffer,
                pe = IgvEvidence.pe_files,
                sr = IgvEvidence.sr_files,
                sample_pe_sr = IgvEvidence.updated_sample_pe_sr,
                reference = reference,
                reference_index = reference_index,
                is_snv_indel = is_snv_indel,
                igv_docker = igv_docker,
                runtime_attr_igv = runtime_attr_igv
        }
    }

    call IntegrateIgvPlots{
        input:
            igv_tar = Igv.tar_gz_pe,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_run_igv
    }

    output {
        File tar_gz_pe = IntegrateIgvPlots.plot_tar
    }
}

task GetSamples {
    input {
        File ped_file
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ped_file, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        Array[String] family_list = read_lines("family_list.txt")
    }

    command <<<
        set -euo pipefail

        cut -f1 ~{ped_file} | sort | uniq > family_list.txt
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

task GetSamplePeSr {
    input {
        File ped_file
        File sample_pe_sr
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
        File sample_pe_sr_list = "sample_pe_sr_list.txt"
    }

    command <<<
        set -euo pipefail

        # Get samples from ped file
        cut -f2 ~{ped_file} > ped_samples.txt

        # Filter PE/SR list to only include samples in ped file
        head -n 1 ~{sample_pe_sr} > sample_pe_sr_list.txt
        while read sample; do
            grep "^$sample" ~{sample_pe_sr} >> sample_pe_sr_list.txt || true
        done < ped_samples.txt
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

task GetFamilySample {
    input {
        File ped_file
        String family
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ped_file, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        Array[String] family_sample_list = read_lines("family_samples.txt")
    }

    command <<<
        set -euo pipefail

        awk -v fam="~{family}" '$1==fam {print $2}' ~{ped_file} > family_samples.txt
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

task IntegrateIgvPlots{
    input {
        Array[File] igv_tar
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(igv_tar, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(50 + input_size * 2),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File plot_tar = "~{prefix}_igv_plots.tar.gz"
    }

    command <<<
        set -euo pipefail

        mkdir ~{prefix}_igv_plots
        for tarfile in ~{sep=" " igv_tar}; do
            tar -xzf "$tarfile"
            mv pe_igv_plots/*  ~{prefix}_igv_plots/
        done < ~{write_lines(igv_tar)};
        tar -czf ~{prefix}_igv_plots.tar.gz ~{prefix}_igv_plots
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