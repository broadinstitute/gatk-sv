version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectQcPerSample {
    input {
        File vcf
        String prefix

        File samples_list

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_override_collect_vids_per_sample
        RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists
    }

    String output_prefix = "~{prefix}.per_sample_qc"

    call CollectVidsPerSample {
        input:
            vcf = vcf,
            samples_list = samples_list,
            prefix = output_prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_override_collect_vids_per_sample
    }

    call MergeShardedPerSampleVidLists {
        input:
            tarballs = [CollectVidsPerSample.vid_lists_tarball],
            samples_list = samples_list,
            prefix = output_prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_override_merge_sharded_per_sample_vid_lists
    }

    output {
        File vid_lists = MergeShardedPerSampleVidLists.merged_tarball
    }
}

task CollectVidsPerSample {
    input {
        File vcf
        File samples_list
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String outdirprefix = "~{prefix}_perSample_VIDs"
    Float input_size = size([vcf, samples_list], "GiB")

    command <<<
        set -euo pipefail

        mkdir -p ~{outdirprefix}

        gq_fmt=$(bcftools view -h ~{vcf} | grep -q '##FORMAT=<ID=GQ' && echo '%GQ' || echo 'NA')

        bcftools view -S ~{samples_list} ~{vcf} \
            | bcftools +fill-tags -- -t AC \
            | bcftools view -i 'AC>0' \
            | bcftools query -f "[%SAMPLE\t%ID\t%GT\t${gq_fmt}\n]" \
            | awk '{OFS="\t"; gt=$3; gsub(/\|/,"/",gt); gq=$4; if(gq=="."){gq=99}; split(gt,a,"/"); x_miss=(a[1]=="."); y_miss=(a[2]=="."); if(x_miss&&y_miss){gt="./."}else if(x_miss){if(a[2]+0>0) gt="0/1"; else gt="0/0"}else if(y_miss){if(a[1]+0>0) gt="0/1"; else gt="0/0"}else{x=a[1]+0;y=a[2]+0; if(x==0&&y==0) gt="0/0"; else if(x==0||y==0) gt="0/1"; else gt="1/1"}; print $1,$2,gt,gq}' \
            | awk -v outprefix="~{outdirprefix}" '$3 != "0/0" && $3 != "./." {OFS="\t"; print $2, $3, $4 >> outprefix"/"$1".VIDs_genotypes.txt" }'

        for FILE in ~{outdirprefix}/*.VIDs_genotypes.txt; do
            gzip -f "$FILE"
            rm -f "$FILE"
        done

        cd ~{outdirprefix} && \
            tar -czvf ../~{outdirprefix}.tar.gz *.VIDs_genotypes.txt.gz && \
            cd -
    >>>

    output {
        File vid_lists_tarball = "~{outdirprefix}.tar.gz"
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(10.0 + input_size),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task MergeShardedPerSampleVidLists {
    input {
        Array[File] tarballs
        File samples_list
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tarballs, "GiB")

    command <<<
        set -euo pipefail

        mkdir "~{prefix}_perSample_VID_lists"
        mkdir shards

        while read i tarball_path; do
            mkdir "shards/shard_$i"
            tar -xzvf "$tarball_path" --directory "shards/shard_$i"/
        done < <( awk -v OFS="\t" '{ print NR, $1 }' ~{write_lines(tarballs)} )

        # Verify tarballs produced files
        N_VID_FILES=$(find shards/ -name "*.VIDs_genotypes.txt.gz" | wc -l)
        if [ "$N_VID_FILES" -eq 0 ]; then
            echo "ERROR: No VID files found in extracted tarballs" >&2
            exit 1
        fi
        echo "Found $N_VID_FILES VID files across shards"

        OUTDIR="~{prefix}_perSample_VID_lists"
        export OUTDIR

        merge_sample() {
            local sample="$1"
            find shards/ -name "${sample}.VIDs_genotypes.txt.gz" \
            | xargs -I {} zcat {} \
            | sort -Vk1,1 -k2,2n -k3,3n \
            | gzip -c \
            > "${OUTDIR}/${sample}.VIDs_genotypes.txt.gz"
        }
        export -f merge_sample

        cat ~{samples_list} | xargs -P $(nproc) -n 1 bash -c 'set -euo pipefail; merge_sample "$1"' _

        # Verify output is non-empty for at least one sample
        N_NONEMPTY=$(find "${OUTDIR}" -name "*.VIDs_genotypes.txt.gz" -size +20c | wc -l)
        if [ "$N_NONEMPTY" -eq 0 ]; then
            echo "ERROR: All merged VID files are empty" >&2
            exit 1
        fi
        echo "$N_NONEMPTY non-empty VID files produced"

        tar -czvf "~{prefix}_perSample_VID_lists.tar.gz" "~{prefix}_perSample_VID_lists"
    >>>

    output {
        File merged_tarball = "~{prefix}_perSample_VID_lists.tar.gz"
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(10 + input_size * 5),
        cpu_cores: 2,
        preemptible_tries: 3,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
