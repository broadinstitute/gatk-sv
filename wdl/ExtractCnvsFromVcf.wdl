version 1.0

workflow ExtractCnvsFromVcf {
	input {
		File vcf               # Input VCF
		File vcf_idx           # Input VCF index file

		String gatk_docker     # Docker image path for GATK
	}
	
	call ExtractCnvs {
		input:
			vcf  					= vcf,
			vcf_idx 	    = vcf_idx,
			gatk_docker 	= gatk_docker
	}

	output {
		File extracted_vcf	 		= ExtractCnvs.extracted_vcf
		File extracted_vcf_idx 	= ExtractCnvs.extracted_vcf_idx
	}
}

task ExtractCnvs {
	input {
		File vcf
		File vcf_idx
		String gatk_docker
	}

	command <<<
		set -eu -o pipefail

		python <<CODE
import pysam
import copy

fin = pysam.VariantFile("~{vcf}")
fout = pysam.VariantFile("extracted.vcf", 'w', header=fin.header)

for rec in fin:
	if rec.info.get("SVTYPE") == "CNV":
		rec_del = copy.copy(rec)
		rec_del.id = rec.id + "_DEL"
		rec_del.info["SVTYPE"] = "DEL"
		rec_del.alts = ("<DEL>",)
		
		rec_dup = copy.copy(rec)
		rec_dup.id = rec.id + "_DUP"
		rec_dup.info["SVTYPE"] = "DUP"
		rec_dup.alts = ("<DUP>",)
		
		fout.write(rec_del)
		fout.write(rec_dup)
	else:
		fout.write(rec)

fin.close()
fout.close()
CODE

			bgzip extracted.vcf

			tabix extracted.vcf.gz
    >>>

	output {
		File extracted_vcf = "extracted.vcf.gz"
		File extracted_vcf_idx = "extracted.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 4 HDD"
		docker: gatk_docker
	}
}