version 1.0

import "Structs.wdl"

workflow IntegratePavVariantsFromHaplotypes {
    input {
        File h1_vcf
        File h2_vcf
        File un_vcf
        String sample

        File? monitoring_script

        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_integrate_variants_in_haplotypes
    }

    call IntegrateVariantsFromHaplotyes{
        input:
            h1_vcf = h1_vcf,
            h2_vcf = h2_vcf,
            un_vcf = un_vcf,
            sample = sample,
            monitoring_script = monitoring_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_integrate_variants_in_haplotypes
    }

    output{
        File output_vcf = IntegrateVariantsFromHaplotyes.inte_vcf
        File output_vcf_idx = IntegrateVariantsFromHaplotyes.inte_vcf_idx
    }
}

task IntegrateVariantsFromHaplotyes {
  input {
    File h1_vcf
    File h2_vcf
    File un_vcf
    String sample
    File? monitoring_script  
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<

    python3 <<CODE

    import argparse
    import gzip

    def open_vcf(file_path):
        return gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r")

    def parse_vcf(file_path):
        header = []
        variants = {}
        with open_vcf(file_path) as f:
            for line in f:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        columns = line.strip().split("\t")
                        sample_col = columns[-1]
                    header.append(line)
                else:
                    cols = line.strip().split("\t")
                    key = tuple(cols[:5])  # CHROM, POS, ID, REF, ALT
                    variants[key] = cols
        return header, variants

    def merge_variants(h1_vars, h2_vars, un_vars):
        all_keys = set(h1_vars.keys()) | set(h2_vars.keys()) | set(un_vars.keys())
        merged = {}
        for key in all_keys:
            sources = []
            if key in h1_vars: sources.append("h1")
            if key in h2_vars: sources.append("h2")
            if key in un_vars: sources.append("un")
            if "h1" in sources and "h2" in sources:
                gt = "1|1"
                record = h1_vars[key]
            elif "h1" in sources:
                gt = "1|0"
                record = h1_vars[key]
            elif "h2" in sources:
                gt = "0|1"
                record = h2_vars[key]
            elif "un" in sources:
                gt = "0/1"
                record = un_vars[key]
            else:
                continue
            # Replace GT field (assume FORMAT is like GT:...)
            format_field = record[8]
            old_sample_fields = record[9].split(":")
            new_sample_fields = [gt] + old_sample_fields[1:] if len(old_sample_fields) > 1 else [gt]
            record[9] = ":".join(new_sample_fields)
            merged[key] = record
        return merged

    def write_merged_vcf(header, variants, output_path, new_sample_name):
        with open(output_path, "w") as out:
            for line in header:
                if line.startswith("#CHROM"):
                    parts = line.strip().split("\t")
                    parts[-1] = new_sample_name
                    out.write("\t".join(parts) + "\n")
                else:
                    out.write(line)
            for record in sorted(variants.values(), key=lambda x: (x[0], int(x[1]))):
                out.write("\t".join(record) + "\n")

    parser = argparse.ArgumentParser(description="Merge h1/h2/un VCFs with proper phasing.")
    parser.add_argument("--h1", required=True, help="Haplotype 1 VCF")
    parser.add_argument("--h2", required=True, help="Haplotype 2 VCF")
    parser.add_argument("--un", required=True, help="Unphased VCF")
    parser.add_argument("-s", "--sample", required=True, help="New sample name for output VCF")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")

    header_h1, h1_vars = parse_vcf("~{h1_vcf}")
    _, h2_vars = parse_vcf("~{h2_vcf}")
    _, un_vars = parse_vcf("~{un_vcf}")
    merged = merge_variants(h1_vars, h2_vars, un_vars)
    write_merged_vcf(header_h1, merged, "PAV_~{sample}.vcf.gz", "~{sample}")


    CODE

    tabix -p vcf "PAV_~{sample}.vcf.gz"
  >>>

  output {
    File inte_vcf = "PAV_~{sample}.vcf.gz"
    File inte_vcf_idx = "PAV_~{sample}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(h1_vcf, "GiB")*6),
    disk_gb: 15 + ceil(size(h1_vcf, "GiB")*6),
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

