version 1.0
import "Structs.wdl"
import "ProcessLrVsSrStat.wdl" as ProcessLrVsSrStat



workflow ProcessLrVsSrStatList {
  input {
    Array[File] full_vcf_list
    Array[File] tp_vcf_list
    Array[File] SVID_GC_list

    Array[String] contig_list

    File vcf2bed_py
    File add_GC_R
    File calcu_stat_R

    String output_prefix
    Boolean related = false   # default is false

    String sv_base_mini_docker
  }

  scatter (i in range(length(full_vcf_list))){
    call ProcessLrVsSrStat.ProcessLrVsSrStat as ProcessLrVsSrStat{
      input:
        full_vcf = full_vcf_list[i],
        tp_vcf = tp_vcf_list[i],
        SVID_GC = SVID_GC_list[i],

        vcf2bed_py = vcf2bed_py,
        add_GC_R = add_GC_R,
        calcu_stat_R = calcu_stat_R,

        contig_list = contig_list,
        related = related,

        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call MergeStatTable as merge_full_stat{
    input:
      freq_files = ProcessLrVsSrStat.full_stat, 
      output_name = output_prefix,
      docker_image = sv_base_mini_docker
  }

  call MergeStatTable as merge_tp_stat{
    input:
      freq_files = ProcessLrVsSrStat.tp_stat, 
      output_name = "~{output_prefix}.tp",
      docker_image = sv_base_mini_docker
  }

  call MergeFillAndTpTables as merge_full_and_tp_stat {
    input:
      file_a = merge_full_stat.merged_table,
      file_b = merge_tp_stat.merged_table,
      output_prefix = output_prefix,
      docker_image = sv_base_mini_docker
  }

  output {
    File merged_stat = merge_full_and_tp_stat.merged_table
  }
}



task MergeStatTable {
  input {
    Array[File] freq_files   # list of input files
    String output_prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    Rscript -e '

      out_file <- "~{output_prefix}.stat"

      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr", repos = "http://cran.us.r-project.org")
      }
      library(dplyr)

      tables <- lapply("~{freq_files}", function(f) {
        df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        required_cols <- c("SVTYPE", "SVLEN_bin", "GC", "AF_bin", "Freq")
        if (!all(required_cols %in% colnames(df))) {
          stop(paste("File", f, "is missing required columns."))
        }
        df
      })

      merged <- bind_rows(tables) %>%
        group_by(SVTYPE, SVLEN_bin, GC, AF_bin) %>%
        summarise(Freq = sum(Freq, na.rm = TRUE), .groups = "drop")

      write.table(merged, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    '
  >>>

  output {
    File merged_table = "~{output_prefix}.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5 ,
    disk_gb: 10 ,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task MergeFillAndTpTables {
  input {
    File file_a
    File file_b
    String output_prefix = "merged_two_freq.tsv"
    String docker_image
    RuntimeAttr? runtime_attr_override
    
    }


  command <<<
    Rscript -e '
      file_a <- "~{file_a}"
      file_b <- "~{file_b}"
      out_file <- "~{output_prefix}.stat"

      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr", repos = "http://cran.us.r-project.org")
      }

      library(dplyr)

      # Read input files
      df_a <- read.table(file_a, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      df_b <- read.table(file_b, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

      required_cols <- c("SVTYPE", "SVLEN_bin", "GC", "AF_bin", "Freq")
      if (!all(required_cols %in% colnames(df_a))) {
        stop("File A is missing required columns.")
      }
      if (!all(required_cols %in% colnames(df_b))) {
        stop("File B is missing required columns.")
      }

      # Rename Freq columns
      df_a <- df_a %>% rename(All = Freq)
      df_b <- df_b %>% rename(TP = Freq)

      # Merge by key columns
      merged <- full_join(df_a, df_b,
                          by = c("SVTYPE", "SVLEN_bin", "GC", "AF_bin")) %>%
                arrange(SVTYPE, SVLEN_bin, GC, AF_bin)

      write.table(merged, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    '
  >>>

  output {
    File merged_table = "~{output_prefix}.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5 ,
    disk_gb: 10 ,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
