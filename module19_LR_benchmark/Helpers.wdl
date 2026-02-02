version 1.0

import "Structs.wdl"

task AddFilter {
    input {
        File vcf
        File vcf_idx
        String filter_name
        String filter_description
        String filter_expression
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -h \
            ~{vcf} \
        | grep "^##" > header.txt
        
        echo '##FILTER=<ID=~{filter_name},Description="~{filter_description}">' >> header.txt
        
        bcftools view \
            -h \
            ~{vcf} \
        | grep "^#CHROM" >> header.txt
        
        bcftools reheader \
            -h header.txt \
            ~{vcf} \
        | bcftools filter \
            --mode + \
            -s ~{filter_name} \
            -e '~{filter_expression}' \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File flagged_vcf = "~{prefix}.vcf.gz"
        File flagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AddInfo {
    input {
        File vcf
        File vcf_idx
        String tag_id
        String tag_value
        String tag_description
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=~{tag_id},Number=1,Type=String,Description="~{tag_description}">' > header.lines

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t~{tag_value}\n' \
            ~{vcf} \
        | bgzip -c > annotations.txt.gz
        
        tabix -s1 -b2 -e2 annotations.txt.gz

        bcftools annotate -h header.lines -a annotations.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/~{tag_id} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateVariantAttributes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=allele_length'; then
            echo '##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=allele_type'; then
            echo '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">' >> new_headers.txt
        fi

        bcftools annotate \
            -h new_headers.txt \
            ~{vcf} \
            -Oz -o temp.vcf.gz
        tabix -p vcf temp.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/allele_length\t%INFO/allele_type\n' \
            temp.vcf.gz \
        | awk -F'\t' '{
            split($4, alleles, ",")
            ref_length = length($3)            
            alt_len = length($4)

            calc_length = alt_len - ref_length
            calc_type = "SNV"
            if (alt_len > ref_len) {
                calc_type = "INS"
            } else if (alt_len < ref_len) {
                calc_type = "DEL"
            }

            allele_length = ($5 == ".") ? calc_length : $5
            allele_type = ($6 == ".") ? calc_type : $6

            print $1"\t"$2"\t"$3"\t"$4"\t"allele_length"\t"allele_type
        }' \
            | bgzip -c > annot.txt.gz
        
        tabix -s1 -b2 -e2 annot.txt.gz

        bcftools annotate -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/allele_length,INFO/allele_type \
            temp.vcf.gz \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String allele_type
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) \
            | sed -e "s/#//g" \
            > ~{allele_type}.bed

        bedtools closest \
            -wo \
            -a <(sort -k1,1 -k2,2n ~{bed_a}) \
            -b <(sort -k1,1 -k2,2n ~{bed_b}) \
            >> ~{allele_type}.bed
    >>>

    output {
        File output_bed = "~{allele_type}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_a, "GB") + size(bed_b, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CheckSampleConsistency {
    input {
        Array[File] vcfs
        Array[File] vcfs_idx
        Array[String] sample_ids
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' sample_ids} | sort > requested_samples.txt

        vcfs_array=(~{sep=' ' vcfs})
        
        for vcf in "${vcfs_array[@]}"; do            
            bcftools query -l "$vcf" | sort > vcf_samples.txt
            
            comm -23 requested_samples.txt vcf_samples.txt > missing.txt
            
            if [ -s missing.txt ]; then
                echo "ERROR: The following samples are missing from $vcf:"
                cat missing.txt
                exit 1
            fi
        done
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatAlignedTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
import csv

input_files = "~{sep=',' tsvs}".split(',')
output_filename = "aligned_unsorted.tsv"
header_filename = "~{prefix}.header.txt"
fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]

all_keys = set()
for f in input_files:
    with open(f, 'r') as fh:
        line = fh.readline().strip()
        if not line: continue
        parts = line.split('\t')
        if len(parts) > 5:
            keys = parts[5:]
            all_keys.update(keys)

sorted_keys = sorted(list(all_keys))
master_header = fixed_cols + sorted_keys

with open(header_filename, 'w') as hout:
    for k in sorted_keys:
        hout.write(k + "\n")

with open(output_filename, 'w') as out:
    out.write("\t".join(master_header) + "\n")
    
    for f in input_files:
        with open(f, 'r') as fh:
            header_line = fh.readline().strip()
            if not header_line: 
                continue
            
            file_cols = header_line.split('\t')
            col_map = {name: i for i, name in enumerate(file_cols)}
            
            for line in fh:
                parts = line.strip().split('\t')
                if not parts: continue
                
                out_row = []
                for target_col in master_header:
                    if target_col in col_map:
                        try:
                            val = parts[col_map[target_col]]
                            out_row.append(val)
                        except IndexError:
                            out_row.append(".")
                    else:
                        out_row.append(".")
                
                out.write("\t".join(out_row) + "\n")
CODE
    
        tail -n +2 aligned_unsorted.tsv | sort -k1,1 -k2,2n > ~{prefix}.tsv
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
        File merged_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        Boolean preserve_header = false
        Boolean skip_sort = false
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ "~{preserve_header}" == "true" ]; then
            head -n 1 ~{tsvs[0]} > ~{prefix}.tsv
            tail -n +2 -q ~{sep=' ' tsvs} >> ~{prefix}.tsv
        elif [ "~{skip_sort}" == "true" ]; then
            cat ~{sep=' ' tsvs} > ~{prefix}.tsv
        else
            cat ~{sep=' ' tsvs} | sort -k1,1 -k2,2n > ~{prefix}.tsv
        fi
    >>>

    output {
        File concatenated_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcfs_idx
        Boolean merge_sort = true
        String prefix = "concat"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String merge_flag = if merge_sort then "--allow-overlaps" else ""

    command <<<
            set -euo pipefail
            
            VCFS_FILE="~{write_lines(vcfs)}"

            bcftools concat \
                ~{merge_flag} \
                --file-list ${VCFS_FILE} \
                -Oz -o "~{prefix}.vcf.gz"
            
            tabix -p vcf -f "~{prefix}.vcf.gz"
    >>>

    output {
            File concat_vcf = "~{prefix}.vcf.gz"
            File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatVcfsLR {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Boolean merge_sort = false
        Boolean remove_dup = true
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    String merge_flag = if merge_sort then "--allow-overlaps" else ""
    Float input_size = size(vcfs, "GB")
    Float compression_factor = 10.0
    Float base_disk_gb = 20.0
    Float base_mem_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + compression_factor * input_size,
        disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            -a ~{merge_flag} \
            --file-list ~{write_lines(vcfs)} \
            -Oz -o merged.tmp.vcf.gz

        tabix -p vcf merged.tmp.vcf.gz

        if [[ ~{remove_dup} == "true" ]]; then
            bcftools norm \
                -d exact \
                -Oz -o ~{prefix}.vcf.gz \
                merged.tmp.vcf.gz
        else
            mv merged.tmp.vcf.gz ~{prefix}.vcf.gz
        fi

        tabix -p vcf ~{prefix}.vcf.gz

    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx =  "~{prefix}.vcf.gz.tbi"
    }
}

task ConvertToSymbolic {
    input {
        File vcf
        File vcf_idx
        String prefix
        Boolean drop_genotypes = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -f '%INFO/allele_type\n' ~{vcf} | sort -u > allele_types.txt

        if [ "~{drop_genotypes}" == "true" ]; then
            bcftools view \
                -G \
                ~{vcf} \
            | python3 /opt/gnomad-lr/scripts/helpers/symbalts.py \
                --input - \
                --output ~{prefix}.vcf.gz \
                --types allele_types.txt
        else
            python3 /opt/gnomad-lr/scripts/helpers/symbalts.py \
                --input ~{vcf} \
                --output ~{prefix}.vcf.gz \
                --types allele_types.txt
        fi
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task DropVcfFields {
    input {
        File vcf
        File vcf_idx
        String drop_fields
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x ~{drop_fields} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File dropped_vcf = "~{prefix}.vcf.gz"
        File dropped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractSample {
    input {
        File vcf
        File vcf_idx
        String sample
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -s ~{sample} \
            --min-ac 1 \
            ~{default="" extra_args} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}

task ExtractVcfAnnotations {
    input {
        File vcf
        File vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        Boolean add_header_row = false
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import sys

vcf = pysam.VariantFile("~{vcf}")
orig = pysam.VariantFile("~{original_vcf}")

new_keys = sorted(list(set(vcf.header.info.keys()) - set(orig.header.info.keys())))

with open("~{prefix}.header.txt", "w") as out:
    for k in new_keys:
        out.write(k + "\n")

with open("~{prefix}.annotations.tsv", "w") as out:
    
    if "~{add_header_row}" == "true":
        fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]
        full_header = fixed_cols + new_keys
        out.write("\t".join(full_header) + "\n")

    for record in vcf:
        alts = ",".join(record.alts) if record.alts else "."
        rid = record.id if record.id else "."
        row = [record.chrom, str(record.pos), record.ref, alts, rid]

        for k in new_keys:
            if k in record.info:
                val = record.info[k]
                if isinstance(val, bool):
                    row.append("1" if val else "0")
                elif isinstance(val, (list, tuple)):
                    row.append(",".join(map(str, val)))
                else:
                    row.append(str(val))
            else:
                row.append(".")
        
        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
        File annotations_header = "~{prefix}.header.txt"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(original_vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FillGenotypesFromUnphased {
    input {
        File phased_vcf
        File phased_vcf_idx
        File unphased_vcf
        File unphased_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import sys

phased_in = pysam.VariantFile("~{phased_vcf}")
unphased_in = pysam.VariantFile("~{unphased_vcf}")

if 'SOURCE' in unphased_in.header.info and 'SOURCE' not in phased_in.header.info:
    phased_in.header.info.add('SOURCE', number=unphased_in.header.info['SOURCE'].number, type=unphased_in.header.info['SOURCE'].type, description=unphased_in.header.info['SOURCE'].description)

for f in unphased_in.header.filters:
    if f not in phased_in.header.filters:
        phased_in.header.filters.add(f, unphased_in.header.filters[f].description)

out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=phased_in.header)

for record in phased_in:
    match = None
    for cand in unphased_in.fetch(record.chrom, record.start, record.stop):
        if cand.id == record.id:
            match = cand
            break
    
    if match is None:
        sys.stderr.write(f"Mismatch: Variant {record.chrom}:{record.pos}-{record.ref}-{record.alts} not found in unphased VCF.\n")
        out.write(record)
        continue

    if 'SOURCE' in match.info:
        record.info['SOURCE'] = match.info['SOURCE']

    if list(match.filter):
        record.filter.clear()
        for f in match.filter:
            record.filter.add(f)
    
    for sample in record.samples:
        u_gt = match.samples[sample]['GT']
        p_gt = record.samples[sample]['GT']
        if u_gt == (0, 0):
            is_target = False
            if p_gt == (None, None):
                is_target = True
            elif len(p_gt) == 2:
                if (p_gt[0] is None and p_gt[1] == 0) or (p_gt[0] == 0 and p_gt[1] is None):
                    is_target = True
            if is_target:
                record.samples[sample]['GT'] = (0, 0)
    
    out.write(record)

out.close()
CODE
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filled_vcf = "~{prefix}.vcf.gz"
        File filled_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(phased_vcf, "GB") + size(unphased_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FinalizeToDir {
    input {
        Array[File] files
        String outdir
        File? keyfile
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    command <<<
        set -euo pipefail

        cat ~{write_lines(files)} | gsutil -m cp -I "~{gcs_output_dir}"
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(files, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FinalizeToFile {
    input {
        File file
        String outdir
        String? name
        File? keyfile
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(file, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GetHailMTSize {
    input {
        String mt_uri
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        tot_size=$(gsutil -m du -sh ~{mt_uri} | awk -F '    ' '{ print $1 }')

        python3 <<CODE > mt_size.txt
import sys

size = "$tot_size".split()[0]
unit = "$tot_size".split()[1]

def convert_to_gib(size, unit):
    size_dict = {"KiB": 2**10, "MiB": 2**20, "GB": 2**30, "TiB": 2**40}
    return float(size) * size_dict[unit] / size_dict["GB"]

size_in_gib = convert_to_gib(size, unit)
print(size_in_gib)
CODE
    >>>

    output {
        Float mt_size = read_lines('mt_size.txt')[0]
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: 25,
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeHeaderLines {
    input {
        Array[File] header_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=' ' header_files} \
            | sort -u > ~{prefix}.txt
    >>>

    output {
        File merged_header = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(header_files, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String? contig
        String? extra_args
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools merge \
            -Oz -o ~{prefix}.vcf.gz \
            ~{if defined(contig) then "-r " + contig else ""} \
            ~{if defined(extra_args) then extra_args else ""} \
            -l ~{write_lines(vcfs)}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}

task RenameVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String id_format = "%CHROM-%POS-%REF-%ALT"
        Boolean strip_chr = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools annotate --set-id '~{id_format}' ~{vcf} -Oz -o temp_renamed.vcf.gz
        
        if [ "~{strip_chr}" == "true" ]; then
            bcftools view -h temp_renamed.vcf.gz > header.txt
            bcftools view -H temp_renamed.vcf.gz \
                | awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $3); print}' \
                | cat header.txt - \
                | bgzip -c > ~{prefix}.vcf.gz
            rm temp_renamed.vcf.gz header.txt
        else
            mv temp_renamed.vcf.gz ~{prefix}.vcf.gz
        fi
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: 8 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10, 
        preemptible_tries: 2, 
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ResetVcfFilters {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x FILTER ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reset_vcf = "~{prefix}.vcf.gz"
        File reset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task RevertSymbolicAlleles {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/helpers/revert_symbalts.py \
            --annotated ~{annotated_vcf} \
            --original ~{original_vcf} \
            --output ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.vcf.gz"
        File reverted_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotated_vcf, "GB") + size(original_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SplitMultiallelics {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import sys

vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("temp.vcf", "w", header=vcf_in.header)

for record in vcf_in:
    if len(record.alts) <= 1:
        vcf_out.write(record)
        continue
    
    ids = record.id.split(';')
    for i, alt_seq in enumerate(record.alts):
        parts = ids[i].split('_')        
        new_rec = record.copy()
        new_rec.chrom = parts[0]
        new_rec.pos = int(parts[1])
        new_rec.ref = parts[2]
        new_rec.alts = (parts[3],)
        new_rec.id = ids[i]
        target_allele_idx = i + 1
        
        for sample in record.samples:
            old_gt = record.samples[sample]['GT']
            new_gt = []
            for allele in old_gt:
                if allele is None:
                    new_gt.append(None)
                elif allele == target_allele_idx:
                    new_gt.append(1)
                else:
                    new_gt.append(0)
            new_rec.samples[sample]['GT'] = tuple(new_gt)
        
        vcf_out.write(new_rec)

vcf_in.close()
vcf_out.close()
CODE

        bcftools sort temp.vcf -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File split_vcf = "~{prefix}.vcf.gz"
        File split_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SplitVcfIntoShards {
    input {
        File input_vcf
        File input_vcf_idx
        Int variants_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir chunks
        bcftools view -h ~{input_vcf} > chunks/header.vcf
        bcftools view -H ~{input_vcf} | split -l ~{variants_per_shard} - chunks/body_

        for body in chunks/body_*; do
        chunk_name=chunks/~{prefix}_$(basename "$body")
        cat chunks/header.vcf "$body" | bgzip -c > "${chunk_name}.vcf.gz"
        tabix -p vcf -f "${chunk_name}.vcf.gz"
        done
    
    >>>

    output {
        Array[File] split_vcfs = glob("chunks/*.vcf.gz")
        Array[File] split_vcf_idxes = glob("chunks/*.vcf.gz.tbi")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(input_vcf,"GB")*2.5) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task StripGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -G \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File stripped_vcf = "~{prefix}.vcf.gz"
        File stripped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetTsvToContig {
    input {
        File tsv
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -v contig="~{contig}" '$1 == contig' ~{tsv} > ~{prefix}.tsv
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfByArgs {
    input {
        File vcf
        File vcf_idx
        String? include_args
        String? exclude_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            ~{if defined(include_args) then "-i '~{include_args}'" else ""} \
            ~{if defined(exclude_args) && !defined(exclude_args) then "-e '~{exclude_args}'" else ""} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfByLength {
    input {
        File vcf
        File vcf_idx
        String length_field = "allele_length"
        String? locus
        Int? min_length
        Int? max_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String size_filter = if defined(min_length) && defined(max_length) then 'abs(INFO/~{length_field})>=~{min_length} && abs(INFO/~{length_field})<=~{max_length}' else if defined(min_length) then 'abs(INFO/~{length_field})>=~{min_length}' else if defined(max_length) then 'abs(INFO/~{length_field})<=~{max_length}' else '1==1'

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include "~{size_filter}" \
            ~{if defined(locus) then "--regions ~{locus}" else ""} \
            -Oz -o ~{prefix}.vcf.gz
                
        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToCalled {
    input {
        File vcf
        File? vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if ~{!defined(vcf_idx)}; then
            tabix -p vcf ~{vcf}
        fi

        bcftools view \
            -i 'ALT != "."' \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToContig {
    input {
        File vcf
        File? vcf_idx
        String contig
        String? args_string
        Boolean drop_genotypes = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        if ~{!defined(vcf_idx)}; then
            tabix -p vcf ~{vcf}
        fi

        bcftools view \
            -r ~{contig} \
            ~{if drop_genotypes then "-G" else ""} \
            ~{if defined(args_string) then args_string else ""} \
            ~{vcf} \
            -Oz -o ~{prefix}.~{contig}.vcf.gz
        tabix -p vcf -f ~{prefix}.~{contig}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.~{contig}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToSampleList {
    input {
        File vcf
        File vcf_idx
        Array[String] samples
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > samples.txt <<EOF
~{sep='\n' samples}
EOF

        bcftools view \
            --samples-file samples.txt \
            --min-ac 1 \
            --regions ~{contig} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SwapSampleIds {
    input {
        File vcf
        File vcf_idx
        File sample_swap_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > current_samples.txt

        awk 'FNR==NR {swap[$1]=$2; next} {if ($1 in swap) print swap[$1]; else print $1}' \
            ~{sample_swap_list} current_samples.txt > new_samples.txt

        bcftools reheader \
            --samples new_samples.txt \
            ~{vcf} \
        > ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File swapped_vcf = "~{prefix}.vcf.gz"
        File swapped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
