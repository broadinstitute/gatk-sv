version 1.0

import "Structs.wdl"

workflow GenotypeGenomicDisorderRegions {
  input {
    Array[String] batch_names
    Array[File] rd_files
    Array[File] median_files
    Array[File] depth_sepcutoff_files

    Array[File] vcfs
    File ped_file
    File genomic_disorder_regions_bed
    File par_bed

    String sv_pipeline_docker
    RuntimeAttr? runtime_generate_median_geno
    RuntimeAttr? runtime_revise_vcf
  }
  scatter (i in range(length(batch_names))) {
    call RunRdTest {
      input:
        batch_name = batch_names[i],
        rdtest_bed = genomic_disorder_regions_bed,
        rd_file = rd_files[i],
        rd_index = rd_files[i] + ".tbi",
        median_file = median_files[i],
        depth_sepcutoff = depth_sepcutoff_files[i],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_generate_median_geno
    }
    call ReviseGenomicDisorderRegions {
      input:
        batch_name = batch_names[i],
        rdtest_tar = RunRdTest.out,
        vcf = vcfs[i],
        ped_file = ped_file,
        genomic_disorder_regions_bed = genomic_disorder_regions_bed,
        par_bed = par_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_revise_vcf
    }
  }
  output{
    Array[File] rdtest_out = RunRdTest.out

    Array[File] revised_records_vcf = ReviseGenomicDisorderRegions.revised_records_vcf
    Array[File] revised_records_index = ReviseGenomicDisorderRegions.revised_records_index
    Array[File] original_records_vcf = ReviseGenomicDisorderRegions.original_records_vcf
    Array[File] original_records_index = ReviseGenomicDisorderRegions.original_records_index
    Array[File] subtracted_vcf = ReviseGenomicDisorderRegions.subtracted_vcf
    Array[File] subtracted_index = ReviseGenomicDisorderRegions.subtracted_index
  }
}

task RunRdTest {
  input{
    String batch_name
    File rdtest_bed
    File rd_file
    File rd_index
    File median_file
    File depth_sepcutoff
    Int large_size_cutoff = 1000000
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 15,
                               disk_gb: ceil(40.0 + size(rd_file, "GiB") * 4),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    # Inject one sample from the batch into the 5th column
    SAMPLE=$(awk -F'\t' '{ if (NR==1) {print $1} }' ~{median_file})
    awk -F'\t' -v OFS='\t' -v s="$SAMPLE" '{print $1,$2,$3,$4,s,$5}' ~{rdtest_bed} > intervals.bed
    mkdir rdtest_~{batch_name}/
    Rscript /opt/RdTest/RdTest.R \
      -v TRUE -g TRUE -p TRUE \
      -r ~{depth_sepcutoff} \
      -b intervals.bed \
      -c ~{rd_file} \
      -m ~{median_file} \
      -n ~{batch_name} \
      -s ~{large_size_cutoff} -o rdtest_~{batch_name}
    tar czvf rdtest_~{batch_name}.tar.gz rdtest_~{batch_name}/
  >>>
  output{
    File out = "rdtest_~{batch_name}.tar.gz"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ReviseGenomicDisorderRegions {
  input{
    String batch_name
    File rdtest_tar
    File vcf
    File ped_file
    File genomic_disorder_regions_bed
    File par_bed
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    tar xzf ~{rdtest_tar}
    RDTEST_DIR=$(basename ~{rdtest_tar} .tar.gz)
    OUTBASE="~{batch_name}.revise_gdr"
    python ~{default="/opt/src/sv-pipeline/scripts/revise_genomic_disorder_regions.py" script} \
      --vcf ~{vcf} \
      --median-geno $RDTEST_DIR/~{batch_name}.median_geno \
      --ped-file ~{ped_file} \
      --region-bed ~{genomic_disorder_regions_bed} \
      --par-bed ~{par_bed} \
      --out $OUTBASE
    mkdir tmp
    bcftools sort -T tmp/ $OUTBASE.new_revised_records.unsorted.vcf.gz -Oz -o $OUTBASE.new_revised_records.vcf.gz
    tabix $OUTBASE.new_revised_records.vcf.gz
    tabix $OUTBASE.original_revised_records.vcf.gz
    tabix $OUTBASE.subtracted.vcf.gz
  >>>
  output{
    File revised_records_vcf = "~{batch_name}.revise_gdr.new_revised_records.vcf.gz"
    File revised_records_index = "~{batch_name}.revise_gdr.new_revised_records.vcf.gz.tbi"
    File original_records_vcf = "~{batch_name}.revise_gdr.original_revised_records.vcf.gz"
    File original_records_index = "~{batch_name}.revise_gdr.original_revised_records.vcf.gz.tbi"
    File subtracted_vcf = "~{batch_name}.revise_gdr.subtracted.vcf.gz"
    File subtracted_index = "~{batch_name}.revise_gdr.subtracted.vcf.gz.tbi"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
