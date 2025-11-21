version 1.0

workflow DefragmentCnvsVcf {
	input {
		File vcf              								# Input VCF
		File vcf_idx          								# Input VCF index file
		Float? max_dist        								# Maximum distance for merging CNVs
		String prefix													# Prefix for output files

		String sv_pipeline_docker							# Docker image path
	}
	
	call Vcf2Bed {
		input:
			vcf = vcf,
			vcf_idx = vcf_idx,
			prefix = prefix,
			sv_pipeline_docker = sv_pipeline_docker
	}

	call DefragmentCnvs {
		input:
			bed_in = Vcf2Bed.bed_out,
			max_dist = max_dist,
			prefix = prefix,
			sv_pipeline_docker = sv_pipeline_docker
	}

	call DefragmentVcfFromBed {
		input:
			vcf = vcf,
			vcf_idx = vcf_idx,
			bed = DefragmentCnvs.bed_out,
			prefix = prefix,
			sv_pipeline_docker = sv_pipeline_docker
	}

	output {
		File defragmented_cnvs_vcf = DefragmentVcfFromBed.out
		File defragmented_cnvs_vcf_idx = DefragmentVcfFromBed.out_idx
	}
}

task Vcf2Bed {
	input {
		File vcf
		File vcf_idx
		String prefix
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		svtk vcf2bed --info ALL ~{vcf} ~{prefix}.input.bed
	>>>

	output {
		File bed_out = "~{prefix}.input.bed"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 4 HDD"
		docker: sv_pipeline_docker
	}
}

task DefragmentCnvs {
	input {
		File bed_in
		Float? max_dist
		String prefix
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		/opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
			--max-dist ~{if defined(max_dist) then max_dist else "0.25"} \
			~{bed_in} \
			~{prefix}.defrag.bed
	>>>

	output {
		File bed_out = "~{prefix}.defrag.bed"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 4 HDD"
		docker: sv_pipeline_docker
	}
}

task DefragmentVcfFromBed {
	input {
		File vcf
		File vcf_idx
		File bed
		String prefix
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		python <<CODE
import pysam

fin = pysam.VariantFile("~{vcf}")
fout = pysam.VariantFile("~{prefix}.defrag.vcf", 'w', header=fin.header)

# Read the BED file and store the intervals
intervals = []
with open("~{bed}") as bed_file:
		for line in bed_file:
				if line.strip() == "":
						continue
				fields = line.strip().split()
				chrom = fields[0]
				start = int(fields[1])
				end = int(fields[2])
				intervals.append((chrom, start, end))

# Merge genotypes for a sample across multiple records
def merge_genotypes(records, sample):
		merged = None
		for rec in records:
				gt = rec.samples[sample].get("GT")
				if gt is None:
						continue
					
				# Convert missing alleles to -1
				alleles = [allele if allele is not None else -1 for allele in gt]

				# If any record is homozygous alt, choose 1/1
				if all(a == 1 for a in alleles):
						return (1, 1)
				
				# If any record is heterozygous, mark as 0/1
				if 0 in alleles and 1 in alleles:
						merged = (0, 1)
				
				# Else, return ./.
		return merged if merged is not None else (None, None)

# For each defragmented interval, fetch overlapping records and merge them
for chrom, start, end in intervals:
		try:
				records = list(fin.fetch(chrom, start, end))
		except ValueError:
				continue
		if not records:
				continue

		# Determine new boundaries
		new_start = min(rec.pos for rec in records) # Alternative: start
		new_end = max(rec.stop for rec in records)  # Alternative: end

		# Create a new record
		new_rec = fout.new_record()
		new_rec.chrom = chrom
		new_rec.pos = new_start
		new_rec.stop = new_end

		# Populate fields based on first record
		new_rec.id = records[0].id + "_defrag"
		new_rec.ref = records[0].ref
		new_rec.alts = records[0].alts
		new_rec.info.update(records[0].info)

		# Merge genotype information across samples
		for sample in fin.header.samples:
				new_rec.samples[sample]["GT"] = merge_genotypes(records, sample)
		
		fout.write(new_rec)

fin.close()
fout.close()
CODE

		# Compress the output VCF
		bgzip ~{prefix}.defrag.vcf

		# Index the compressed VCF
		tabix ~{prefix}.defrag.vcf.gz
	>>>

	output {
		File out = "~{prefix}.defrag.vcf.gz"
		File out_idx = "~{prefix}.defrag.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 4 HDD"
		docker: sv_pipeline_docker
	}
}
