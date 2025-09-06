version 1.0

import "Structs.wdl"

workflow Igv {
    input {
        String family
        Array[String] samples
        File cram_crai_list
        File varfile
        Boolean requester_pays
        Int igv_max_window
        Boolean file_localization
        File reference
        File reference_index
        Boolean is_snv_indel
        String buffer
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_localize_reads
    }

    if(file_localization){
        call RunIGV_whole_genome_localize {
            input:
                family = family,
                samples = samples,
                cram_crai_list = cram_crai_list,
                varfile = varfile,
                requester_pays = requester_pays,
                igv_max_window = igv_max_window,
                buffer = buffer,
                reference = reference,
                reference_index = reference_index,
                is_snv_indel = is_snv_indel,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    if(!file_localization){
        call RunIGV_whole_genome_parse{
            input:
                family = family,
                samples = samples,
                cram_crai_list = cram_crai_list,
                varfile = varfile,
                igv_max_window = igv_max_window,
                buffer = buffer,
                reference = reference,
                reference_index = reference_index,
                is_snv_indel = is_snv_indel,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    output {
        File? tar_gz_pe = select_first([RunIGV_whole_genome_localize.pe_plots, RunIGV_whole_genome_parse.pe_plots])
    }
}

task RunIGV_whole_genome_localize{
    input {
        String family
        Array[String] samples
        File cram_crai_list
        File varfile
        Boolean requester_pays
        Int igv_max_window
        String buffer
        File reference
        File reference_index
        Boolean is_snv_indel
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([cram_crai_list, varfile], "GB")
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
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
    }

    command <<<
        set -euo pipefail

        mkdir pe_igv_plots

        # Process variant file
        cat ~{varfile} | awk 'NR>1 {print $1":"$2"-"$3"," $4"," $5}' > variant_info.csv

        # Create CRAM files list
        while IFS=$'\t' read -r sample crai cram; do
            echo "$sample $cram $crai" >> crams.txt
        done < ~{cram_crai_list}

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
            
            # Check if window is too large
            if [[ $((window_end - window_start)) -gt ~{igv_max_window} ]]; then
                # Split into two windows
                mid=$((start + (end - start) / 2))
                echo "${chrom}:${window_start}-$((mid + ~{buffer}))" >> regions_${i}.txt
                echo "${chrom}:$((mid - ~{buffer}))-${window_end}" >> regions_${i}.txt
            else
                echo "${region}" >> regions_${i}.txt
            fi

            # Create IGV script
            cat > pe.${i}.txt << EOF
new
genome hg38
load ~{reference}
EOF

            # Add CRAM files for samples in this family
            for sample in ~{sep=" " samples}; do
                cram_file=$(awk -v s="$sample" '$1==s {print $2}' crams.txt)
                if [[ -n "$cram_file" ]]; then
                    echo "load $cram_file" >> pe.${i}.txt
                fi
            done

            # Add regions and snapshots
            while read region; do
                echo "goto $region" >> pe.${i}.txt
                echo "collapse" >> pe.${i}.txt
                echo "sort position" >> pe.${i}.txt
                echo "snapshot pe_igv_plots/${variant_id}_${family}.png" >> pe.${i}.txt
            done < regions_${i}.txt

            echo "exit" >> pe.${i}.txt

            # Run IGV (placeholder - would need actual IGV setup)
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

task RunIGV_whole_genome_parse{
    input {
        String family
        Array[String] samples
        File cram_crai_list
        File varfile
        Int igv_max_window
        String buffer
        File reference
        File reference_index
        Boolean is_snv_indel
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([cram_crai_list, varfile], "GB")
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
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
    }

    command <<<
        set -euo pipefail

        mkdir pe_igv_plots

        # Process variant file
        cat ~{varfile} | awk 'NR>1 {print $1":"$2"-"$3"," $4"," $5}' > variant_info.csv

        # Create CRAM files list
        while IFS=$'\t' read -r sample crai cram; do
            echo "$sample $cram $crai" >> crams.txt
        done < ~{cram_crai_list}

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
            
            # Split large regions
            if [[ $((window_end - window_start)) -gt ~{igv_max_window} ]]; then
                echo "${chrom}:${window_start}-$((start + ~{buffer}))" >> regions_${i}.txt
                echo "${chrom}:$((end - ~{buffer}))-${window_end}" >> regions_${i}.txt
            else
                echo "${region}" >> regions_${i}.txt
            fi

            # Create IGV script
            cat > pe.${i}.txt << EOF
new
genome hg38
load ~{reference}
EOF

            # Add CRAM files for samples in this family
            for sample in ~{sep=" " samples}; do
                cram_file=$(awk -v s="$sample" '$1==s {print $2}' crams.txt)
                if [[ -n "$cram_file" ]]; then
                    echo "load $cram_file" >> pe.${i}.txt
                fi
            done

            # Add regions and snapshots
            while read region; do
                echo "goto $region" >> pe.${i}.txt
                echo "collapse" >> pe.${i}.txt
                echo "sort position" >> pe.${i}.txt
                echo "snapshot pe_igv_plots/${variant_id}_${family}.png" >> pe.${i}.txt
            done < regions_${i}.txt

            echo "exit" >> pe.${i}.txt

            # Run IGV (placeholder - would need actual IGV setup)
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