version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow HailMerge {
  input {
    Array[File] vcfs
    String prefix
    String? gcs_project  # REQUIRED
    Boolean? reset_cnv_gts
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_fix_header
  }

  # Concatenate vcfs naively to prevent ClassTooLargeException in Hail
  if (length(vcfs) > 1) {
    call MiniTasks.ConcatVcfs as Preconcat {
      input:
        vcfs=vcfs,
        naive=true,
        generate_index=false,
        outfile_prefix="~{prefix}.preconcat",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_preconcat
    }
  }

  call HailMergeTask {
    input:
      vcfs = [select_first([Preconcat.concat_vcf, vcfs[0]])],
      prefix = prefix,
      gcs_project = select_first([gcs_project]),
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      runtime_attr_override=runtime_override_hail_merge
  }

  call FixHeader {
    input:
      merged_vcf = HailMergeTask.merged_vcf,
      example_vcf = vcfs[0],
      prefix = prefix + ".reheadered",
      reset_cnv_gts = select_first([reset_cnv_gts, false]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override=runtime_override_fix_header
  }

  output {
    File merged_vcf = FixHeader.out
    File merged_vcf_index = FixHeader.out_index
  }
}

task HailMergeTask {
  input {
    Array[File] vcfs
    String prefix
    String gcs_project
    String region = "us-central1"
    String sv_pipeline_hail_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcfs: {
      localization_optional: true
    }
  }

  String cluster_name_prefix="gatk-sv-cluster-"

  RuntimeAttr runtime_default = object {
                                  mem_gb: 6.5,
                                  disk_gb: 100,
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " SSD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_hail_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    cp ~{write_lines(vcfs)} "files.list"

    python <<CODE
import hail as hl
import os
import uuid
from google.cloud import dataproc_v1 as dataproc

cluster_name = "gatk-sv-hail-{}".format(uuid.uuid4())
script_path = "/opt/sv-pipeline/scripts/hailmerge.py"

try:
  print(os.popen("hailctl dataproc start --num-workers 4 --region {} --project {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=1440m {}".format("~{region}", "~{gcs_project}", cluster_name)).read())

  cluster_client = dataproc.ClusterControllerClient(
        client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
  )

  for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
    if cluster.cluster_name == cluster_name:
      cluster_staging_bucket = cluster.config.temp_bucket
      os.popen("gcloud dataproc jobs submit pyspark {} --cluster={} --project {} --files=files.list --region={} --driver-log-levels root=WARN -- {} {}".format(script_path, cluster_name, "~{gcs_project}", "~{region}", cluster_staging_bucket, cluster_name)).read()
      os.popen("gsutil cp -r gs://{}/{}/merged.vcf.bgz .".format(cluster_staging_bucket, cluster_name)).read()
      break

except Exception as e:
  print(e)
  raise
finally:
  os.popen("gcloud dataproc clusters delete --project {} --region {} {}".format("~{gcs_project}", "~{region}", cluster_name)).read()
CODE

  mv merged.vcf.bgz ~{prefix}.vcf.gz
  tabix ~{prefix}.vcf.gz
  >>>

  output {
    File merged_vcf = "~{prefix}.vcf.gz"
    File merged_vcf_index = "~{prefix}.vcf.gz.tbi"
  }
}

task FixHeader {
  input {
    File merged_vcf
    File example_vcf
    String prefix
    Boolean reset_cnv_gts
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + size(merged_vcf, "GB") * 2 + size(example_vcf, "GB")),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " SSD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # Reset to original header
    bcftools view --no-version -h ~{merged_vcf}  | grep -v ^#CHROM > header
    bcftools view --no-version -h ~{example_vcf} | grep -e "^##source" -e "^##ALT" -e "^##CPX_TYPE" >> header
    bcftools view --no-version -h ~{merged_vcf}  | grep ^#CHROM >> header
    bcftools reheader -h header ~{merged_vcf} \
      ~{if reset_cnv_gts then "| gunzip | python /opt/sv-pipeline/04_variant_resolution/scripts/reset_cnv_gts.py stdin stdout | bgzip" else ""} \
      > ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
