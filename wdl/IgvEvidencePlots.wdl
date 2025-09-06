version 1.0

import "Structs.wdl"

workflow Igv {
    input {
        String family
        Array[String] samples
        File varfile
        String buffer
        File pe
        File sr
        File sample_pe_sr
        File reference
        File reference_index
        Boolean is_snv_indel
        String igv_docker
        RuntimeAttr? runtime_attr_igv
    }

    call RunIGV_whole_genome {
        input:
            family = family,
            samples = samples,
            varfile = varfile,
            buffer = buffer,
            pe = pe,
            sr = sr,
            sample_pe_sr = sample_pe_sr,
            reference = reference,
            reference_index = reference_index,
            is_snv_indel = is_snv_indel,
            igv_docker = igv_docker,
            runtime_attr_override = runtime_attr_igv
    }

    output {
        File tar_gz_pe = RunIGV_whole_genome.pe_plots
    }
}

task RunIGV_whole_genome {
    input {
        String family
        Array[String] samples
        File varfile
        String buffer
        File pe
        File sr
        File sample_pe_sr
        File reference
        File reference_index
        Boolean is_snv_indel
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([varfile, pe, sr], "GB")
    Float base_mem_gb = 7.5

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(50 + input_size * 5),
        cpu: 1,
        preemptible: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File pe_plots = "~{family}_pe_igv_plots.tar.gz"
    }

    command <<<
        set -euo pipefail

        mkdir pe_igv_plots

        # Process variant file
        cat ~{varfile} | awk 'NR>1 {print $1":"$2"-"$3"," $4"," $5}' > variant_info.csv

        # Generate IGV scripts for each variant
        i=0
        while IFS=',' read -r region variant_id svtype; do
            i=$((i+1))
            chrom=$(echo "$region" | cut -d':' -f1)
            coords=$(echo "$region" | cut -d':' -f2)
            start=$(echo "$coords" | cut -d'-' -f1)
            end=$(echo "$coords" | cut -d'-' -f2)
            
            # Calculate window with buffer
            window_start=$((start - ~{buffer}))
            window_end=$((end + ~{buffer}))
            
            # Create IGV script
            cat > pe.${i}.txt << EOF
new
genome hg38
load ~{reference}
load ~{pe}
load ~{sr}
goto ${chrom}:${window_start}-${window_end}
collapse
sort position
snapshot pe_igv_plots/${variant_id}_${family}.png
exit
EOF

            # Placeholder for IGV execution
            echo "Would run IGV with script pe.${i}.txt" > pe_igv_plots/${variant_id}_${family}.log

        done < variant_info.csv

        tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
} 