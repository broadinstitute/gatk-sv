import gzip

# Input and output files
in_vcf = "/Users/kjaising/Downloads/cluster_batch.dragen.fmt.vcf.gz"
out_vcf = "/Users/kjaising/Downloads/cluster_batch.dragen.gq.fmt.vcf.gz"

# Default GQ value to use
DEFAULT_GQ = "999"

with gzip.open(in_vcf, 'rt') as f_in, gzip.open(out_vcf, 'wt') as f_out:
    for line in f_in:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                f_out.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            f_out.write(line)
        else:
            fields = line.strip().split('\t')
            format_idx = 8
            format_fields = fields[format_idx].split(':')

            if 'GQ' not in format_fields:
                format_fields.insert(1, 'GQ')
                fields[format_idx] = ':'.join(format_fields)

                for i in range(9, len(fields)):
                    sample_values = fields[i].split(':')
                    sample_values.insert(1, DEFAULT_GQ)
                    fields[i] = ':'.join(sample_values)

            f_out.write('\t'.join(fields) + '\n')

print(f"Created new VCF with GQ fields at: {out_vcf}")
