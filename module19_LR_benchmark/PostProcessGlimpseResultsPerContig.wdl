version 1.0

import "Structs.wdl"
import "MergeVcfs.wdl" as MergeVcfs
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow PostProcessGlimpseResultsPerContig {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_tbi

        File genotype_output_vcf
        File genotype_output_vcf_tbi

        File panel_biallelic_vcf
        File panel_biallelic_vcf_idx

        File convert_to_biallelic_script
        File? monitoring_script

        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_split_ref_panel
        RuntimeAttr? runtime_attr_add_id_to_info_column
        RuntimeAttr? runtime_attr_convert_bubbles_to_biallelic
        RuntimeAttr? runtime_attr_deduplicate_vcf
    }


    call SplitRefPanel{
        input:
            ref_panel_vcf = ref_panel_vcf,
            ref_panel_vcf_tbi = ref_panel_vcf_tbi,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_split_ref_panel
    }

    call AddIdToInfoColumn {
        input:
            a_vcf = genotype_output_vcf,
            a_vcf_idx = genotype_output_vcf_tbi,
            b_vcf = SplitRefPanel.split_id_vcf,
            b_vcf_idx = SplitRefPanel.split_id_vcf_idx,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_add_id_to_info_column
        }

    call LongReadGenotypeTasks.ConvertBubblesToBiallelic{
        input:
            input_vcf = AddIdToInfoColumn.with_id_vcf,
            input_vcf_idx = AddIdToInfoColumn.with_id_vcf_idx,

            panel_biallelic_vcf = panel_biallelic_vcf,
            panel_biallelic_vcf_idx = panel_biallelic_vcf_idx,

            sort = true,
            convert_to_biallelic_script = convert_to_biallelic_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_convert_bubbles_to_biallelic
    }

    call DeDuplicateBiallelicVcfs{
        input:
            biallelic_vcf = ConvertBubblesToBiallelic.biallelic_vcf,
            biallelic_vcf_idx = ConvertBubblesToBiallelic.biallelic_vcf_idx, 
            docker_image = sv_pipeline_base_docker, 
            runtime_attr_override = runtime_attr_deduplicate_vcf
  }


    output{
        File output_biallelic_vcf = DeDuplicateBiallelicVcfs.deduplicated_vcf
        File output_biallelic_vcf_idx = DeDuplicateBiallelicVcfs.deduplicated_vcf_idx
    }
}


#task1: split the multi-allelic sites in the reference panel to bi-allelic, split ID in the info column accordingly
task SplitRefPanel {
  input {
    File ref_panel_vcf
    File ref_panel_vcf_tbi  
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(ref_panel_vcf,'.vcf.gz')

  command <<<

    python3 <<CODE
    import os
    import sys
    import pysam

    def split_vcf(input_vcf, output_vcf):
        vcf_in = pysam.VariantFile(input_vcf, 'r')
        vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

        for record in vcf_in:
            if len(record.alts) > 1:
                id_list = record.info['ID']
                for idx, alt in enumerate(record.alts):
                    new_record = record.copy()
                    # Update ALT and keep only the one allele
                    new_record.alts = (alt,)
                    new_record.info['ID'] = (id_list[idx])
                    # Clear the ID field for the split records (optional)
                    new_record.id = id_list[idx]
                    # Adjust genotypes (if needed) — skipped for simplicity
                    vcf_out.write(new_record)
            else:
                # Just pass through records with a single ALT allele
                record.id =record.info['ID'][0]
                vcf_out.write(record)

        vcf_in.close()
        vcf_out.close()

    split_vcf("~{ref_panel_vcf}", "~{prefix}.split_id.vcf.gz")

    CODE

    tabix -p vcf "~{prefix}.split_id.vcf.gz"
  >>>

  output {
    File split_id_vcf = "~{prefix}.split_id.vcf.gz"
    File split_id_vcf_idx = "~{prefix}.split_id.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 25,
    disk_gb: 50 + ceil(size(ref_panel_vcf, "GiB")*2),
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

#task2: add the ID info in the info column to the output vcf from GLIMPSE
task AddIdToInfoColumn {
  input {
    File a_vcf
    File a_vcf_idx
    File b_vcf
    File b_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(a_vcf,'.vcf.gz')

  command <<<

    python3 <<CODE

    import pysam
    import sys

    def variant_key(record):
        """Create a unique key for matching: (CHROM, POS, REF, ALT)"""
        return (record.chrom, record.pos, record.ref, tuple(record.alts))

    def load_vcf_b_ids(vcf_b_path):
        match_dict = {}
        with pysam.VariantFile(vcf_b_path) as vcf_b:
            for record in vcf_b:
                for alt in record.alts:
                    key = (record.chrom, record.pos, record.ref, alt)
                    match_dict[key] = record.info['ID']
        return match_dict

    def add_ids_from_b_to_a(vcf_a_path, vcf_b_path, output_path):
        # Load match dictionary from VCF B
        match_ids = load_vcf_b_ids(vcf_b_path)
        with pysam.VariantFile(vcf_a_path) as vcf_a:
            # Add new INFO field for MatchID if not already present
            header = vcf_a.header
            if "ID" not in header.info:
                header.info.add("ID", "A", "String", "Variant IDs per ALT allele.")
            with pysam.VariantFile(output_path, 'w', header=header) as vcf_out:
                for record in vcf_a:
                    new_record = record.copy()
                    match_found = False
                    for alt in record.alts:
                        key = (record.chrom, record.pos, record.ref, alt)
                        if key in match_ids:
                            new_record.info['ID'] = match_ids[key]
                            new_record.id = match_ids[key][0]
                        else: 
                            print(key)
                    vcf_out.write(new_record)


    add_ids_from_b_to_a("~{a_vcf}", "~{b_vcf}", "~{prefix}.with_id.vcf.gz")

    CODE

    tabix -p vcf "~{prefix}.with_id.vcf.gz"
  >>>

  output {
    File with_id_vcf = "~{prefix}.with_id.vcf.gz"
    File with_id_vcf_idx = "~{prefix}.with_id.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 30,
    disk_gb: 50 + ceil(size(a_vcf, "GiB") *10),
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

task DeDuplicateBiallelicVcfs{
  input {
    File biallelic_vcf
    File biallelic_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(biallelic_vcf,'.vcf.gz')

  command <<<

    bcftools view -H ~{biallelic_vcf} | cut -f3 | sort | uniq -c | awk '{if ($1>1) print $2}' > duplicated_SVID.tsv

    python3 <<CODE
    #!/usr/bin/env python3

    import argparse
    import pysam
    from collections import defaultdict, Counter

    def parse_svid_list(svid_file):
        with open(svid_file) as f:
            return set(line.strip() for line in f if line.strip())

    def is_alt_gt(gt):
        return gt in {(0, 1), (1, 0), (1, 1)}

    def is_alt_hom_gt(gt):
        return gt in {(1, 1)}

    def select_genotype(gts):
        count = Counter(gts)
        most_common = count.most_common()
        top_count = most_common[0][1]
        top_gts = [gt for gt, c in most_common if c == top_count]
        
        if len(top_gts) == 1:
            return top_gts[0]
        else:
            # Prefer genotype with alt allele if tie
            for gt in top_gts:
                if is_alt_hom_gt(gt):
                    return gt
                elif is_alt_gt(gt):
                    return gt
            return top_gts[0]

    def collapse_vcf_by_svid(input_vcf, svid_file, output_vcf, output_2_vcf):
        svids = parse_svid_list(svid_file)
        vcf_in = pysam.VariantFile(input_vcf)
        vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
        vcf_2_out = pysam.VariantFile(output_2_vcf, 'w', header=vcf_in.header)

        # Group records by ID (SVID)
        grouped = defaultdict(list)
        for rec in vcf_in.fetch():
            if rec.id in svids:
                grouped[rec.id].append(rec)
            else:
                vcf_2_out.write(rec)

        for svid, records in grouped.items():
            new_rec = records[0]  # Take the first as representative
            for sample in vcf_in.header.samples:
                gts = [rec.samples[sample]["GT"] for rec in records]
                chosen_gt = select_genotype(gts)
                new_rec.samples[sample]["GT"] = chosen_gt

            vcf_out.write(new_rec)

        vcf_out.close()
        vcf_2_out.close()
        print(f"Output written to {output_vcf}")



    input_vcf = 'chr22.with_id_biallelic.vcf.gz'
    svid_file = 'duplicated_SVID.tsv'
    output_vcf = 'chr22.with_id_biallelic.dedup.vcf.gz'
    output2_vcf = 'chr22.with_id_biallelic.uniq.vcf.gz'

    collapse_vcf_by_svid("~{biallelic_vcf}", "duplicated_SVID.tsv", "~{prefix}.dedup.vcf.gz", "~{prefix}.uniq.vcf.gz")

    CODE

    tabix -p vcf "~{prefix}.dedup.vcf.gz"
    tabix -p vcf "~{prefix}.uniq.vcf.gz"
    bcftools concat -a -Oz -o "~{prefix}.dedup.sort.vcf.gz"  "~{prefix}.dedup.vcf.gz" "~{prefix}.uniq.vcf.gz"
    tabix -p vcf "~{prefix}.dedup.sort.vcf.gz"

  >>>

  output {
    File deduplicated_vcf =  "~{prefix}.dedup.sort.vcf.gz"
    File deduplicated_vcf_idx = "~{prefix}.dedup.sort.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 30,
    disk_gb: 50 + ceil(size(biallelic_vcf, "GiB") *5),
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


