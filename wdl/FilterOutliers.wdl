version 1.0

# Contact Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"

# Workflow to identify & filter outliers from VCFs after module 03 (random forest)
workflow FilterOutlierSamples {
  input {
    String batch
    Array[File?] vcfs
    Array[String] samples
    Array[String] algorithms
    Int N_IQR_cutoff
    File? outlier_cutoff_table
    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_exclude_outliers
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
  }

  Int num_algorithms = length(algorithms)

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      if (defined(outlier_cutoff_table)) {
        call IdentifyOutliersByCutoffTable {
          input:
            vcf = select_first([vcfs[i]]),
            outlier_cutoff_table = select_first([outlier_cutoff_table]),
            outfile = "${algorithms[i]}_outliers.txt",
            algorithm = algorithms[i],
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_identify_outliers
        }
      }
      call IdentifyOutliersByIQR {
        input:
          vcf = select_first([vcfs[i]]),
          N_IQR_cutoff = N_IQR_cutoff,
          outfile = "${algorithms[i]}_outliers.txt",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_identify_outliers
      }
    }
  }

  # Merge list of outliers
  call CatOutliers {
    input:
      outliers = flatten([select_all(IdentifyOutliersByIQR.outliers_list),select_all(IdentifyOutliersByCutoffTable.outliers_list)]),
      batch = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_outliers
  }

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      call ExcludeOutliers {
        input:
          vcf = select_first([vcfs[i]]),
          outliers_list = CatOutliers.outliers_list,
          outfile = "${batch}.${algorithms[i]}.outliers_removed.vcf.gz",
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_exclude_outliers
      }
    }
  }

  # Write new list of samples without outliers
  call FilterSampleList {
    input:
      original_samples = samples,
      outlier_samples = CatOutliers.outliers_list,
      batch = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_filter_samples
  }

  output {
    Array[File?] vcfs_noOutliers = ExcludeOutliers.vcf_no_outliers
    Array[String] filtered_batch_samples_list = FilterSampleList.filtered_samples_list
    File filtered_batch_samples_file = FilterSampleList.filtered_samples_file
    Array[String] outlier_samples_excluded = CatOutliers.outliers_list
    File outlier_samples_excluded_file = CatOutliers.outliers_file
  }
}

task IdentifyOutliersByIQR {
  input {
    File vcf
    Int N_IQR_cutoff
    String outfile
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File outliers_list = "${outfile}"
  }
  command <<<

    set -euo pipefail
    # Count sv per class per sample
    svtk count-svtypes ~{vcf} svcounts.txt

    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/03_variant_filtering/scripts/get_outliers_from_svcounts.helper.R \
      svcounts.txt \
      ~{N_IQR_cutoff} \
      ~{outfile}
  
  >>>
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

task IdentifyOutliersByCutoffTable {
  input {
    File vcf
    File outlier_cutoff_table
    String outfile
    String algorithm
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File outliers_list = "${outfile}"
  }
  command <<<

    set -euo pipefail
    # Count sv per class per sample
    svtk count-svtypes ~{vcf} svcounts.txt

    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/03_variant_filtering/scripts/get_outliers_from_svcounts.helper_V2.R \
      svcounts.txt \
      ~{outlier_cutoff_table} \
      ~{outfile} \
      ~{algorithm}

  >>>
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

# Exclude outliers from VCF
task ExcludeOutliers {
  input {
    File vcf
    Array[String] outliers_list
    String outfile
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf_no_outliers = "${outfile}"
  }
  command <<<

    set -eu
    OUTLIERS=~{write_lines(outliers_list)}
    if [ $( wc -c < $OUTLIERS ) -gt 1 ]; then
      zcat ~{vcf} | fgrep "#" | fgrep -v "##" \
       | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
       | fgrep -wf $OUTLIERS | cut -f2 \
       > indexes_to_exclude.txt
      zcat ~{vcf} | \
       cut --complement -f$( cat indexes_to_exclude.txt | paste -s -d, ) \
       | vcftools --mac 1 --vcf - --recode --recode-INFO-all --stdout \
       | bgzip -c > ~{outfile}
    else
      cp ~{vcf} ~{outfile}
    fi
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

# Write new list of samples per batch after outlier filtering
task FilterSampleList {
  input {
    Array[String] original_samples
    Array[String] outlier_samples
    String batch
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[String] filtered_samples_list = read_lines("${batch}.post03_outliers_excluded.samples.list")
    File filtered_samples_file = "${batch}.post03_outliers_excluded.samples.list"
  }
  command <<<

    fgrep -wvf ~{write_lines(outlier_samples)} ~{write_lines(original_samples)} > ~{batch}.post03_outliers_excluded.samples.list
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

# Merge outlier sample lists across algorithms
task CatOutliers {
  input {
    Array[File] outliers
    String batch
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[String] outliers_list = read_lines("${batch}.post03_outliers.samples.list")
    File outliers_file = "${batch}.post03_outliers.samples.list"
  }
  command <<<

    set -euo pipefail
    while read file; do
      [ -e "$file" ] || continue
      cat $file
    done < ~{write_lines(outliers)} | sort | uniq > ~{batch}.post03_outliers.samples.list
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
