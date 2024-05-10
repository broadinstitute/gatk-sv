version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as format

workflow FormatVcfForGatkScattered {
  input {
    Array[File] vcfs
  }

  scatter ( vcf in vcfs) {
    call format.FormatVcfForGatk {
      input:
        vcf=vcf
    }
  }
  output {
    Array[File] gatk_formatted_vcf = FormatVcfForGatk.gatk_formatted_vcf
    Array[File] gatk_formatted_vcf_index = FormatVcfForGatk.gatk_formatted_vcf_index
  }
}
