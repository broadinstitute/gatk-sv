version 1.0

import "Structs.wdl"

workflow LiftoverVCFs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File chain_file
    String liftover_ref_version

    String liftover_docker
    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_vcf_to_bed
    RuntimeAttr? runtime_attr_liftover
    RuntimeAttr? runtime_attr_update_vcf
    RuntimeAttr? runtime_attr_sort_index

  }

  scatter(i in range(length(vcfs))) {
    call VCFToBED {
      input:
        vcf = vcfs[i],
        vcf_idx = vcf_idxs[i],
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_vcf_to_bed
    }

    call Liftover {
      input:
        bed = VCFToBED.bed,
        chain_file = chain_file,
        liftover_ref_version = liftover_ref_version,
        docker_image = liftover_docker,
        runtime_attr_override = runtime_attr_liftover
    }

    call UpdateVCF {
      input:
        vcf = vcfs[i],
        vcf_idx = vcf_idxs[i],
        bed = Liftover.lifted_bed,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_update_vcf
    }

    call SortAndIndex {
      input:
        vcf = UpdateVCF.updated_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_sort_index
    }
  }

    output {
        Array[File] lifted_sorted_vcfs = SortAndIndex.sorted_vcf
        Array[File] lifted_sorted_tbis = SortAndIndex.sorted_tbi
    }
}



task VCFToBED {
    input {
        File vcf
        File vcf_idx
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        python3 <<CODE

        #!/usr/bin/env python3
        import sys
        import gzip

        def open_file(filename):
            if filename.endswith(".gz"):
                return gzip.open(filename, "rt")
            else:
                return open(filename, "r")

        def vcf_to_bed(vcf_file, bed_file):
            with open_file(vcf_file) as fin, open(bed_file, "w") as fout:
                fout.write("#chr\tstart\tend\tid\tref\talt\n")
                for line in fin:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 5:
                        continue
                    chrom, pos, vid, ref, alts = fields[:5]
                    try:
                        start = int(pos)
                        end = start + len(ref) - 1
                    except ValueError:
                        continue
                    fout.write(f"{chrom}\t{start}\t{end}\t{vid}\t{ref}\t{alts}\n")

        vcf_to_bed("~{vcf}", "~{prefix}.bed")

        CODE
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 15 + ceil(size(vcf, "GiB") *5),
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
}

task Liftover {
    input {
        File bed
        File chain_file
        String liftover_ref_version
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bed, ".bed")

    command <<<
        set -euo pipefail
        liftOver ~{bed} ~{chain_file} ~{prefix}.~{liftover_ref_version}.liftover.bed ~{prefix}.~{liftover_ref_version}.liftover.unmap
    >>>

    output {
        File lifted_bed = "~{prefix}.~{liftover_ref_version}.liftover.bed"
        File unmap = "~{prefix}.~{liftover_ref_version}.liftover.unmap"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 15 + ceil(size(bed, "GiB") *5),
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
}

task UpdateVCF {
    input {
        File vcf
        File vcf_idx
        File bed
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bed, ".bed")
    command <<<
        set -euo pipefail

        python3 <<CODE

        #!/usr/bin/env python3
        import sys
        import gzip

        def open_file(filename, mode="rt"):
            """Open plain text or gzipped files."""
            if filename.endswith(".gz"):
                return gzip.open(filename, mode)
            return open(filename, mode)

        def load_bed(bed_file):
            """Load BED into dict {ID: (chr, pos)}"""
            bed_dict = {}
            with open_file(bed_file, "rt") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 4:
                        continue
                    chrom, pos, _, vid = cols[:4]
                    bed_dict[vid] = (chrom, pos)
            return bed_dict

        def update_vcf(vcf_file, bed_dict, output_file):
            """Update CHROM and POS based on BED coordinates."""
            with open_file(vcf_file, "rt") as fin, gzip.open(output_file, "wt") as fout:
                for line in fin:
                    if line.startswith("#"):
                        fout.write(line)
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 8:
                        continue
                    vid = cols[2]
                    if vid in bed_dict:
                        new_chr, new_pos = bed_dict[vid]
                        cols[0] = new_chr
                        cols[1] = new_pos
                    fout.write("\t".join(cols) + "\n")

        vcf_file = "~{vcf}"
        bed_file = "~{bed}"
        output_file = "~{prefix}.vcf.gz"

        bed_dict = load_bed(bed_file)
        print(f"Loaded {len(bed_dict)} entries from BED file.")
        update_vcf(vcf_file, bed_dict, output_file)
        print(f"Updated VCF written to: {output_file}")

        CODE

    >>>

    output {
        File updated_vcf = "~{prefix}.vcf.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 15 + ceil(size(vcf, "GiB") *5),
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
}

task SortAndIndex {
    input {
        File vcf
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix    = basename(vcf, ".vcf.gz")
    command <<<
        set -euo pipefail

        bcftools sort ~{vcf} -Oz -o ~{prefix}.sorted.vcf.gz
        tabix -p vcf ~{prefix}.sorted.vcf.gz
    >>>

    output {
        File sorted_vcf = "~{prefix}.sorted.vcf.gz"
        File sorted_tbi = "~{prefix}.sorted.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 15 + ceil(size(vcf, "GiB") *5),
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
}
