version 1.0

import "SVConcordance.wdl" as conc
import "Utils.wdl" as util
import "TasksMakeCohortVcf.wdl" as mini_tasks

workflow SVConcordancePerSample {
  input {
    File cohort_eval_vcf
    File cohort_truth_vcf
    File sample_ids_file
    String prefix
    String sv_base_mini_docker
  }

  Array[String] sample_ids = read_lines(sample_ids_file)

  scatter (sample_id in sample_ids) {
    call util.SubsetVcfBySamplesList as SubsetEval {
      input:
        vcf = cohort_eval_vcf,
        vcf_idx = cohort_eval_vcf + ".tbi",
        list_of_samples = write_lines([sample_id]),
        sv_base_mini_docker = sv_base_mini_docker
    }

    call util.SubsetVcfBySamplesList as SubsetTruth {
      input:
        vcf = cohort_truth_vcf,
        vcf_idx = cohort_truth_vcf + ".tbi",
        list_of_samples = write_lines([sample_id]),
        sv_base_mini_docker = sv_base_mini_docker
    }

    call conc.SVConcordance {
      input:
        eval_vcf = SubsetEval.vcf_subset,
        truth_vcf = SubsetTruth.vcf_subset,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call utils.TarFiles {
    input:
      files=flatten([SVConcordance.concordance_vcf, SVConcordance.concordance_vcf_index]),
      prefix="~{prefix}.sample_concordance_vcfs"
  }

  output {
    File persample_concordance_vcf_tar = TarFiles.out
  }
}
