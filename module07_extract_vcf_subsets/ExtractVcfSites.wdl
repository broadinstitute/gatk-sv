## Copyright Broad Institute, 2022
## 
##
## Extract VCF subsets by sample list or extract sites only.
## Optionally shard VCF into smaller chunks by variant count.
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"

workflow ExtractVcfSites {
    input {
        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[String] contigs
        File? sample_list
        Int? shard_counts
        String sv_pipeline_docker
        String sv_base_mini_docker
        String sv_pipeline_hail_docker
        RuntimeAttr? runtime_attr_override_extract_subset_vcf
    }

    scatter(i in range(length(vcf_list))) {
        # If no sharding, process VCF as whole
        if (!defined(shard_counts)) {
            call ExtractSubsetSamples {
                input:
                    vcf = vcf_list[i],
                    vcf_idx = vcf_idx_list[i],
                    sample_list = sample_list,
                    midfix = contigs[i],
                    extract_sites_only = !defined(sample_list),
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_override_extract_subset_vcf
            }
        }

        # If sharding, first shard the VCF then extract samples/sites from each shard
        if (defined(shard_counts)) {
            call ShardVcf {
                input:
                    vcf = vcf_list[i],
                    vcf_idx = vcf_idx_list[i],
                    shard_counts = select_first([shard_counts]),
                    midfix = contigs[i],
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_override_extract_subset_vcf
            }

            scatter(j in range(length(ShardVcf.shard_vcfs))) {
                call ExtractSubsetSamples as ExtractFromShard {
                    input:
                        vcf = ShardVcf.shard_vcfs[j],
                        vcf_idx = ShardVcf.shard_vcf_idxs[j],
                        sample_list = sample_list,
                        midfix = contigs[i] + ".shard_~{j}",
                        extract_sites_only = !defined(sample_list),
                        sv_pipeline_docker = sv_pipeline_docker,
                        runtime_attr_override = runtime_attr_override_extract_subset_vcf
                }
            }
        }
    }

    Array[File] subset_vcf_list_unsharded = select_all(ExtractSubsetSamples.out_vcf)
    Array[Array[File]] subset_vcf_list_sharded_nested = select_all(ExtractFromShard.out_vcf)
    Array[File] subset_vcf_list_sharded = flatten(subset_vcf_list_sharded_nested)
    Array[File] subset_vcf_list_all = flatten([subset_vcf_list_unsharded, subset_vcf_list_sharded])

    Array[File] subset_vcf_idx_list_unsharded = select_all(ExtractSubsetSamples.out_vcf_idx)
    Array[Array[File]] subset_vcf_idx_list_sharded_nested = select_all(ExtractFromShard.out_vcf_idx)
    Array[File] subset_vcf_idx_list_sharded = flatten(subset_vcf_idx_list_sharded_nested)
    Array[File] subset_vcf_idx_list_all = flatten([subset_vcf_idx_list_unsharded, subset_vcf_idx_list_sharded])

    output {
        Array[File] subset_vcf_list = subset_vcf_list_all
        Array[File] subset_vcf_idx_list = subset_vcf_idx_list_all
    }
}



task ShardVcf {
    input {
        File vcf
        File vcf_idx
        Int shard_counts
        String midfix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 20.0
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + (input_size * 3.0)),
        cpu_cores: 2,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(vcf, '.vcf.gz')
    command <<<
        set -eu -o pipefail

        python3 << 'PYTHON_EOF'
import gzip
import subprocess
import sys

vcf_file = "~{vcf}"
shard_size = ~{shard_counts}
prefix = "~{prefix}.~{midfix}"

shard_num = 0
header_lines = []
variant_count = 0
current_shard_gz = None

print(f"Sharding {vcf_file} with {shard_size} variants per shard", file=sys.stderr)

with gzip.open(vcf_file, 'rt') as fin:
    for line in fin:
        if line.startswith('#'):
            header_lines.append(line)
        else:
            # Start new shard if needed
            if variant_count % shard_size == 0:
                if current_shard_gz is not None:
                    current_shard_gz.close()
                    # Index the previous shard
                    prev_shard = f"{prefix}.shard_{shard_num}.vcf.gz"
                    subprocess.run(['tabix', '-p', 'vcf', prev_shard], check=True)
                
                shard_num += 1
                current_shard_file = f"{prefix}.shard_{shard_num}.vcf.gz"
                current_shard_gz = gzip.open(current_shard_file, 'wt')
                # Write header to new shard
                for h in header_lines:
                    current_shard_gz.write(h)

            current_shard_gz.write(line)
            variant_count += 1

# Close and index the last shard
if current_shard_gz is not None:
    current_shard_gz.close()
    last_shard = f"{prefix}.shard_{shard_num}.vcf.gz"
    subprocess.run(['tabix', '-p', 'vcf', last_shard], check=True)

print(f"Created {shard_num} shards from {variant_count} total variants", file=sys.stderr)

PYTHON_EOF

    >>>

    output {
        Array[File] shard_vcfs = glob("~{prefix}.shard_*.vcf.gz")
        Array[File] shard_vcf_idxs = glob("~{prefix}.shard_*.vcf.gz.tbi")
    }

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task ExtractSubsetSamples {
    input {
        File vcf
        File vcf_idx
        File? sample_list
        String midfix
        Boolean extract_sites_only
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: ceil(base_disk_gb + (input_size * 2.0)),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(vcf, '.vcf.gz')
    command <<<
        set -eu -o pipefail

        if ~{extract_sites_only}; then
            # Extract sites only (no genotypes)
            bcftools view -G ~{vcf} \
            | bgzip > ~{prefix}.~{midfix}.vcf.gz
        else
            # Extract samples if sample_list is provided, otherwise keep all genotypes
            if [[ -f ~{sample_list} ]]; then
                bcftools view -S ~{sample_list} ~{vcf} \
                | bgzip > ~{prefix}.~{midfix}.vcf.gz
            else
                bcftools view ~{vcf} \
                | bgzip > ~{prefix}.~{midfix}.vcf.gz
            fi
        fi

        tabix -p vcf ~{prefix}.~{midfix}.vcf.gz

    >>>

    output {
        File out_vcf = "~{prefix}.~{midfix}.vcf.gz"
        File out_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
    }

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

