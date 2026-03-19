version 1.0

import "Structs.wdl"

workflow SplitVCFToHaplotypesWorkflow {
  input {
    File input_vcf
    File input_vcf_idx
    String output_prefix
    String output_midfix
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_split_vcf_to_haplotypes
  }

  call SplitVCFToHaplotypes {
    input:
      input_vcf = input_vcf,
      input_vcf_idx = input_vcf_idx,
      output_prefix = output_prefix,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_vcf_to_haplotypes
  }

  output {
    File hap_vcf = SplitVCFToHaplotypes.output_vcf
  }
}

task SplitVCFToHaplotypes {
  input {
    File input_vcf
    File input_vcf_idx
    String output_prefix
    String output_type = "z"  # v, z, or b

    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_override

  }

  command <<<
    set -euo pipefail

    cat << 'EOF' > split_vcf_haplotypes.py
    #!/usr/bin/env python3

    import argparse
    import pysam
    import sys
    import os

    def parse_args():
        parser = argparse.ArgumentParser(
            description="Split phased diploid genotypes into haplotypes"
        )
        parser.add_argument("-i", "--input", required=True)
        parser.add_argument("-o", "--output", required=True)
        parser.add_argument("--output-type", choices=["v","z","b"], default="z")
        return parser.parse_args()

    def main():
        args = parse_args()

        mode_map = {"v": "w", "z": "wz", "b": "wb"}
        mode = mode_map[args.output_type]

        vcf_in = pysam.VariantFile(args.input)
        old_samples = list(vcf_in.header.samples)

        new_header = pysam.VariantHeader()

        # Copy header records
        for rec in vcf_in.header.records:
            if rec.key != "SAMPLE":
                new_header.add_record(rec)

        # Add haplotype samples
        for s in old_samples:
            new_header.add_sample(f"{s}_hap1")
            new_header.add_sample(f"{s}_hap2")

        vcf_out = pysam.VariantFile(args.output, mode, header=new_header)

        for rec in vcf_in:
            new_rec = vcf_out.new_record(
                contig=rec.contig,
                start=rec.start,
                stop=rec.stop,
                alleles=rec.alleles,
                id=rec.id,
                qual=rec.qual,
                filter=rec.filter.keys()
            )

            for key in rec.info:
                new_rec.info[key] = rec.info[key]

            for s in old_samples:
                gt = rec.samples[s].get("GT")

                if gt is None or None in gt:
                    h1, h2 = (None, None)
                else:
                    if len(gt) == 2:
                        h1, h2 = gt
                    else:
                        h1, h2 = (None, None)

                new_rec.samples[f"{s}_hap1"]["GT"] = (h1,)
                new_rec.samples[f"{s}_hap2"]["GT"] = (h2,)

            vcf_out.write(new_rec)

        vcf_in.close()
        vcf_out.close()

    if __name__ == "__main__":
        main()
    EOF

    python3 split_vcf_haplotypes.py \
        -i ~{input_vcf} \
        -o ~{output_prefix}.vcf.gz \
        --output-type ~{output_type}
  >>>

  output {
    File output_vcf = "~{output_prefix}.vcf.gz"
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(input_vcf,"GiB"))*2,
        disk_gb: 20 + ceil(size(input_vcf,"GiB"))*2,
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
