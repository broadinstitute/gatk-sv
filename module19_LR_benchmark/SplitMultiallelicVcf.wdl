version 1.0

import "Structs.wdl"

# Split multi-allelic VCF records into individual biallelic records.
# For each split record:
#   - ID is appended with _1, _2, ... for each alt allele
#   - Per-allele INFO fields (Number=A in header, or comma-count matches n_alts)
#     are split to their respective value; TYPE="." values are inferred from
#     relative lengths of REF vs ALT
#   - Non-per-allele INFO fields are preserved as-is
#   - Genotypes are recoded: target allele -> 1, ref -> 0, other alts -> .

workflow SplitMultiallelicVcfWorkflow {
  input {
    File vcf_gz
    File vcf_idx
    String output_prefix
    String docker_image

    RuntimeAttr? runtime_attr_override
  }

  call SplitMultiallelicVcf {
    input:
      vcf_gz          = vcf_gz,
      vcf_idx         = vcf_idx,
      output_prefix   = output_prefix,
      docker_image    = docker_image,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File output_vcf_gz  = SplitMultiallelicVcf.output_vcf_gz
    File output_vcf_idx = SplitMultiallelicVcf.output_vcf_idx
  }
}

task SplitMultiallelicVcf {
  input {
    File   vcf_gz
    File   vcf_idx
    String output_prefix
    String docker_image

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores:        1,
    mem_gb:           4,
    disk_gb:          50,
    boot_disk_gb:     10,
    preemptible_tries: 3,
    max_retries:      1
  }

  RuntimeAttr runtime_attr =
    select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import gzip
import re
import sys

def infer_type(ref, alt):
    """Infer variant type from ref and alt allele lengths."""
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    elif len(ref) == len(alt):
        return "mnp"
    elif len(alt) > len(ref):
        return "ins"
    else:
        return "del"


def recode_gt(gt_str, allele_idx):
    """Recode a GT string for a specific alt allele index (1-based).
    
    - allele == allele_idx  -> 1
    - allele == 0           -> 0
    - allele == other alt   -> . (missing)
    """
    phased = '|' in gt_str
    sep = '|' if phased else '/'
    alleles = gt_str.split(sep)
    new_alleles = []
    for a in alleles:
        if a == '.':
            new_alleles.append('.')
        elif a == '0':
            new_alleles.append('0')
        elif int(a) == allele_idx:
            new_alleles.append('1')
        else:
            new_alleles.append('.')
    return sep.join(new_alleles)


def split_multiallelic(input_vcf, output_vcf):
    open_in = gzip.open if input_vcf.endswith('.gz') else open

    # Collect per-allele INFO field IDs from header (Number=A)
    per_allele_info = set()

    with open_in(input_vcf, 'rt') as fin, open(output_vcf, 'wt') as fout:
        for line in fin:
            # --- Header lines ---
            if line.startswith('##'):
                fout.write(line)
                m = re.match(r'##INFO=<ID=([^,]+),Number=A,', line)
                if m:
                    per_allele_info.add(m.group(1))
                continue

            if line.startswith('#'):
                fout.write(line)
                continue

            # --- Variant lines ---
            fields = line.rstrip('\n').split('\t')
            chrom, pos, vid, ref, alt_str = fields[0], fields[1], fields[2], fields[3], fields[4]
            qual, filt, info_str = fields[5], fields[6], fields[7]

            alt_alleles = alt_str.split(',')
            n_alts = len(alt_alleles)

            has_samples = len(fields) >= 10
            format_str  = fields[8] if has_samples else None
            samples     = fields[9:] if has_samples else []

            # Helper: decide if an INFO key/value pair is per-allele
            def is_per_allele(key, val):
                if key in per_allele_info:
                    return True
                # Heuristic: if value has exactly n_alts comma-separated parts
                if val is not None and len(val.split(',')) == n_alts:
                    return True
                return False

            # Parse INFO into ordered list of (key, value_or_None) pairs
            info_pairs = []
            for item in info_str.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_pairs.append((k, v))
                else:
                    info_pairs.append((item, None))

            # Biallelic records: only fix TYPE="." if needed, then pass through
            if n_alts == 1:
                new_info_items = []
                for k, v in info_pairs:
                    if k == 'TYPE' and v == '.':
                        new_info_items.append(f'TYPE={infer_type(ref, alt_alleles[0])}')
                    elif v is not None:
                        new_info_items.append(f'{k}={v}')
                    else:
                        new_info_items.append(k)
                out_fields = fields[:7] + [';'.join(new_info_items)] + fields[8:]
                fout.write('\t'.join(out_fields) + '\n')
                continue

            # Multi-allelic: emit one record per alt allele
            for i, alt in enumerate(alt_alleles, start=1):
                # Unique ID: original ID + _<allele_index>
                new_id = f"{vid}_{i}"

                # Build per-split INFO
                new_info_items = []
                for k, v in info_pairs:
                    if v is None:
                        new_info_items.append(k)
                    elif is_per_allele(k, v):
                        parts = v.split(',')
                        val_i = parts[i - 1] if (i - 1) < len(parts) else '.'
                        # Infer TYPE when value is "."
                        if k == 'TYPE' and val_i == '.':
                            val_i = infer_type(ref, alt)
                        new_info_items.append(f'{k}={val_i}')
                    else:
                        new_info_items.append(f'{k}={v}')
                new_info = ';'.join(new_info_items)

                # Recode genotypes
                if has_samples:
                    fmt_fields = format_str.split(':')
                    gt_idx = fmt_fields.index('GT') if 'GT' in fmt_fields else None
                    new_samples = []
                    for sample in samples:
                        s_fields = sample.split(':')
                        if gt_idx is not None and gt_idx < len(s_fields):
                            s_fields[gt_idx] = recode_gt(s_fields[gt_idx], i)
                        new_samples.append(':'.join(s_fields))
                    out_fields = [chrom, pos, new_id, ref, alt, qual, filt, new_info,
                                  format_str] + new_samples
                else:
                    out_fields = [chrom, pos, new_id, ref, alt, qual, filt, new_info]

                fout.write('\t'.join(out_fields) + '\n')


split_multiallelic("~{vcf_gz}", "~{output_prefix}.split_multiallelic.vcf")
print("Split complete.", file=sys.stderr)
PYEOF
    bgzip "~{output_prefix}.split_multiallelic.vcf"
    tabix -p vcf "~{output_prefix}.split_multiallelic.vcf.gz"
  >>>

  output {
    File output_vcf_gz  = "~{output_prefix}.split_multiallelic.vcf.gz"
    File output_vcf_idx = "~{output_prefix}.split_multiallelic.vcf.gz.tbi"
  }

  runtime {
    cpu:           select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:        select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks:         "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,    default_attr.boot_disk_gb])
    docker:        docker_image
    preemptible:   select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:    select_first([runtime_attr.max_retries,       default_attr.max_retries])
  }
}
