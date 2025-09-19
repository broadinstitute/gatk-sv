version 1.0
import "Structs.wdl"


workflow GTConcordanceWorkflow {
  input {
    File pav_vcf
    File kg_vcf
    File pg_vcf
    File pav_vcf_tbi
    File kg_vcf_tbi
    File pg_vcf_tbi
    Array[String] contigs
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_concordance
  }


  scatter (contig in contigs) {

    call SplitVCFs as split_input_vcf {
      input:
        input_vcf = pav_vcf,
        input_vcf_tbi = pav_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }


   call SplitVCFs as split_kg_vcf {
      input:
        input_vcf = kg_vcf,
        input_vcf_tbi = kg_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }

   call SplitVCFs as split_pg_vcf {
      input:
        input_vcf = pg_vcf,
        input_vcf_tbi = pg_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }

    call RunConcordance {
      input:
        pav_contig = split_input_vcf.contig_vcf,
        kg_contig = split_kg_vcf.contig_vcf,
        pg_contig = split_pg_vcf.contig_vcf,
        docker_image = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_concordance
    }
  }

  output {
    Array[File] kg_outputs = RunConcordance.kg_output
    Array[File] pg_outputs = RunConcordance.pg_output
  }
}

task SplitVCFs {
  input {
    File input_vcf
    File input_vcf_tbi
    String contig
    String docker_image
  }

  String prefix = basename(input_vcf, ".vcf.gz") 
  
  command <<<
    bcftools view -r ~{contig} -Oz -o ~{prefix}.~{contig}.vcf.gz ~{input_vcf}
   >>>

  output {
    File contig_vcf = "~{prefix}.~{contig}.vcf.gz"

  }

  runtime {
    docker: docker_image
    memory: "4G"
    cpu: 1
  }
}


task RunConcordance {
  input {
    File pav_contig
    File kg_contig
    File pg_contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String kg_output = sub(basename(kg_contig), ".vcf.gz", ".SVID_concor")
  String pg_output = sub(basename(pg_contig), ".vcf.gz", ".SVID_concor")

  command <<<

  set -euxo pipefail

  Rscript -e ' 

  calc_genotype_concordance <- function(df) {
    # Extract genotypes from method 1 (cols 5–13) and method 2 (cols 14–22)
    gt1 <- as.character(unlist(df[1, 5:13]))
    gt2 <- as.character(unlist(df[1, 14:22]))

    # Construct table of 9 rows × 2 columns
    gt_table <- data.frame(method1 = gt1, method2 = gt2, stringsAsFactors = FALSE)
    gt_table[,1] = gsub("[|]", "/" , gt_table[,1])

    # Exclude rows with missing genotypes (".", "./.", NA)
    valid_rows <- !(grepl("\\.", gt_table$method1) | grepl("\\.", gt_table$method2))

    gt_table <- gt_table[valid_rows, ]
    effective_gt = nrow(gt_table)

    #if (nrow(gt_table) == 0) {
    #  return(c(effective_gt,NA,NA,NA))
    #}
    if (nrow(gt_table) == 0) {
      return(list(
        effective_gt = effective_gt,
        concordance_exact = NA,
        concordance_altref = NA,
        concordance_nonref = NA
      ))
    }


    ## 1. Exact concordance
    matches_exact <- gt_table$method1 == gt_table$method2
    concordance_exact <- sum(matches_exact) / nrow(gt_table)

    ## 2. Alt vs Ref concordance
    to_altref <- function(gt) {
      if (gt %in% c("0/0")) return("Ref")
      if (gt %in% c("0/1", "1/0", "1/1")) return("Alt")
      return(NA)
    }

    gt_table$method1_altref <- vapply(gt_table$method1, to_altref, character(1))
    gt_table$method2_altref <- vapply(gt_table$method2, to_altref, character(1))

    valid_altref <- !(is.na(gt_table$method1_altref) | is.na(gt_table$method2_altref))
    matches_altref <- gt_table$method1_altref[valid_altref] == gt_table$method2_altref[valid_altref]
    concordance_altref <- if (sum(valid_altref) > 0) sum(matches_altref) / sum(valid_altref) else NA

    ## 3. Non-ref concordance (only rows where at least one method ≠ 0/0)
    nonref_rows <- (gt_table$method1 != "0/0") | (gt_table$method2 != "0/0")
    gt_nonref <- gt_table[nonref_rows, ]

    if (nrow(gt_nonref) == 0) {
      concordance_nonref <- NA
    } else {
      matches_nonref <- gt_nonref$method1 == gt_nonref$method2
      concordance_nonref <- sum(matches_nonref) / nrow(gt_nonref)
    }

    ## Return results
    #return(c(effective_gt,concordance_exact,concordance_altref,concordance_nonref))
      return(list(
          effective_gt = effective_gt,
          concordance_exact = concordance_exact,
          concordance_altref = concordance_altref,
          concordance_nonref = concordance_nonref
        ))


    }

  apply_concordance <- function(df) {
    # Apply the concordance function to each row
    results <- apply(df, 1, function(row) {
      res <- calc_genotype_concordance(as.data.frame(t(row), stringsAsFactors = FALSE))
      unlist(res)
    })

    # Transpose results and convert to data.frame
    results_df <- as.data.frame(t(results), stringsAsFactors = FALSE)

    # Make sure numeric
    results_df <- data.frame(lapply(results_df, as.numeric))

    # Bind to original table
    df_out <- cbind(df, results_df)
    return(df_out)
  }

  # Wrapper to split, process in chunks, and combine results
  apply_concordance_in_chunks <- function(df, chunk_size = 1000) {
    n <- nrow(df)
    idx <- seq(1, n, by = chunk_size)

    results_list <- lapply(idx, function(i) {
      j <- min(i + chunk_size - 1, n)  # end index
      print(i)
      df_chunk <- df[i:j, , drop = FALSE]
      apply_concordance(df_chunk)
    })

    # Combine all chunks
    df_out <- do.call(rbind, results_list)
    return(df_out)
  }

  # ----------------------------
  # Load data
  # ----------------------------
  d1 <- read.table("~{pav_contig}", header=FALSE)
  kg <- read.table("~{kg_contig}", header=FALSE)
  pg <- read.table("~{pg_contig}", header=FALSE)
  for(i in c(10:ncol(pg))){
    pg[,i] = sapply(pg[,i], function(x){strsplit(as.character(x),":")[[1]][1]})
  }


  # ----------------------------
  # Fix column names
  # ----------------------------
  samples <- c("HG00512","HG00513","HG00514",
               "HG00731","HG00732","HG00733",
               "NA19238","NA19239","NA19240")

  colnames(d1)[10:18] <- samples
  colnames(kg)[10:18] <- samples
  colnames(pg)[10:18] <- samples

  # ----------------------------
  # Merge
  # ----------------------------
  pav_kg <- merge(d1[,c(1,2,4,5,10:18)], kg[,c(1,2,4,5,10:18)], by=c("V1","V2","V4","V5"))
  pav_pg <- merge(d1[,c(1,2,4,5,10:18)], pg[,c(1,2,4,5,10:18)], by=c("V1","V2","V4","V5"))

  # ----------------------------
  # Run concordance
  # ----------------------------
  results_kg <- apply_concordance_in_chunks(pav_kg, chunk_size=1000)
  results_pg <- apply_concordance_in_chunks(pav_pg, chunk_size=1000)

  # ----------------------------
  # Save results
  # ----------------------------
  write.table(results_kg, "~{kg_output}", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  write.table(results_pg, "~{pg_output}", sep="\t", col.names=TRUE, row.names=FALSE)
 
  '

  >>>

  output {
    File kg_output = "~{kg_output}"
    File pg_output = "~{pg_output}"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 1, 
      mem_gb: 15, 
      disk_gb: 30,
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
      docker: sv_fst_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }    
}
