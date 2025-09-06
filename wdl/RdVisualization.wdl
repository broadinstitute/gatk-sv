version 1.0

import "Structs.wdl"

workflow RdTestVisualization {
    input {
        String prefix
        File ped_file
        File? fam_ids
        File batch_medianfile
        File batch_bincov
        File bed
        String sv_pipeline_rdtest_docker
        String variant_interpretation_docker
        File outlier_samples
        File sample_batches
        RuntimeAttr? runtime_attr_rdtest
    }

    call CreateRdTestBed {
        input:
            bed = bed,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker
    }

    call RdTest {
        input:
            bed = CreateRdTestBed.bed_output,
            batch_medianfile = batch_medianfile,
            batch_bincov = batch_bincov,
            prefix = prefix,
            outlier_samples = outlier_samples,
            sample_batches = sample_batches,
            ped_file = ped_file,
            fam_ids = fam_ids,
            sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
            runtime_attr_override = runtime_attr_rdtest
    }

    output {
        File Plots = RdTest.plots
    }
}

task CreateRdTestBed {
    input {
        File bed
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size * 1.5),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File bed_output = "~{prefix}.rdtest.bed"
    }

    command <<<
        set -euo pipefail

        # Convert variant bed to RdTest format: chr start end name sample svtype
        # Skip header lines and process each variant
        grep -v "^#" ~{bed} | awk -F'\t' 'NF >= 6 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $5}' > ~{prefix}.rdtest.bed
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

task RdTest {
    input {
        File bed
        File batch_medianfile
        File batch_bincov
        String prefix
        File outlier_samples
        File sample_batches
        File ped_file
        File? fam_ids
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([bed, batch_medianfile, batch_bincov], "GB")
    Float base_mem_gb = 7.5

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(50 + input_size * 10),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File plots = "~{prefix}_rd_plots.tar.gz"
    }

    command <<<
        set -euo pipefail

        # Create directory for plots
        mkdir ~{prefix}_rd_plots

        # Get list of samples from batch file
        cut -f1 ~{sample_batches} | tail -n+2 > samples.txt

        # Filter out outlier samples if provided
        if [[ -s ~{outlier_samples} ]]; then
            grep -v -f ~{outlier_samples} samples.txt > filtered_samples.txt || cp samples.txt filtered_samples.txt
        else
            cp samples.txt filtered_samples.txt
        fi

        # Run RdTest
        Rscript /opt/RdTest/RdTestV2.R \
            -b ~{bed} \
            -n ~{prefix} \
            -m ~{batch_medianfile} \
            -c ~{batch_bincov} \
            -f ~{ped_file} \
            -p TRUE \
            -w filtered_samples.txt \
            -o ~{prefix}_rd_plots/

        # Create tarball
        tar -czf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_pipeline_rdtest_docker
    }
} 