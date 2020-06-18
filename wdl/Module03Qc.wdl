##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/03_variant_filtering/27/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "MasterVcfQc.wdl" as vcf_qc

workflow Module03Qc {
  input {
    File manta_vcf_noOutliers
    File delly_vcf_noOutliers
    File melt_vcf_noOutliers
    File depth_vcf_noOutliers
    File merged_pesr_vcf
    String batch
    File ped_file
    Array[File]? thousand_genomes_tarballs
    Array[File]? hgsv_tarballs
    Array[File]? asc_tarballs
    File? sanders_2015_tarball
    File? collins_2017_tarball
    File? werling_2018_tarball
  }

  call vcf_qc.MasterVcfQc as delly_qc {
    input:
      vcf=delly_vcf_noOutliers,
      ped_file=ped_file,
      prefix="${batch}.delly_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball
  }

  call vcf_qc.MasterVcfQc as manta_qc {
    input:
      vcf=manta_vcf_noOutliers,
      ped_file=ped_file,
      prefix="${batch}.manta_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball
  }

  call vcf_qc.MasterVcfQc as melt_qc {
    input:
      vcf=melt_vcf_noOutliers,
      ped_file=ped_file,
      prefix="${batch}.melt_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball
  }

  call vcf_qc.MasterVcfQc as pesr_qc {
    input:
      vcf=merged_pesr_vcf,
      ped_file=ped_file,
      prefix="${batch}.pesr_merged_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball
  }

  call vcf_qc.MasterVcfQc as depth_qc {
    input:
      vcf=depth_vcf_noOutliers,
      ped_file=ped_file,
      prefix="${batch}.depth_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball
  }

  output {
    File filtered_delly_vcf_qc = delly_qc.sv_vcf_qc_output
    File filtered_manta_vcf_qc = manta_qc.sv_vcf_qc_output
    File filtered_melt_vcf_qc = melt_qc.sv_vcf_qc_output
    File filtered_pesr_vcf_qc = pesr_qc.sv_vcf_qc_output
    File filtered_depth_vcf_qc = depth_qc.sv_vcf_qc_output
  }

}
