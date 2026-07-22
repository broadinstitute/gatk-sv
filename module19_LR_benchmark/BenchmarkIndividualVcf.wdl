version 1.0

import "Structs.wdl"
import "ExtractFileByIndex.wdl" as ExtractFileByIndex
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "BenchmarkIndividualVcfPerContig.wdl" as BenchmarkIndividualVcfPerContig
import "SplitMultiallelicVcf.wdl" as SplitMultiallelicVcf


workflow BenchmarkIndividualVcf{
    input{
        File ref_vcf
        File ref_vcf_idx
        File? ref_filter_vcf
        File? ref_filter_vcf_idx
        File query_vcf
        File query_vcf_idx  
        File? query_filter_vcf
        File? query_filter_vcf_idx

        Array[String] chromosomes

        File anno_script_bash
        File anno_script_helper_R
        File benchmark_script_bash
        File banchmark_script_helper_R


        File repeat_mask
        File simple_repeats
        File segmental_duplicates

        Array[String] sample_ids
        File ref_dict

        Boolean short_read_benchmark = false

        String? truvari_params
        String sv_base_mini_docker
        String sv_pipeline_base_docker

        Boolean run_split_multiallelic_query = true
        Boolean run_split_multiallelic_ref = true

        RuntimeAttr? runtime_attr_split_multiallelic
    }

    String query_prefix = basename(query_vcf, ".vcf.gz")
    String ref_prefix = basename(ref_vcf, ".vcf.gz")

    scatter (index in range(length(chromosomes))){

        call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query{
            input:
                input_vcf = query_vcf,
                input_vcf_index = query_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{query_prefix}.~{chromosomes[index]}.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }

        call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref{
            input:
                input_vcf = ref_vcf,
                input_vcf_index = ref_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }

        if (run_split_multiallelic_query) {
          call SplitMultiallelicVcf.SplitMultiallelicVcf as split_multiallelic_query{
              input:
                  vcf_gz       = extract_chrom_variants_query.chr_vcf,
                  vcf_idx      = extract_chrom_variants_query.chr_vcf_idx,
                  output_prefix = "~{query_prefix}.~{chromosomes[index]}.split",
                  docker_image = sv_pipeline_base_docker,
                  runtime_attr_override = runtime_attr_split_multiallelic
          }
        }

        if (run_split_multiallelic_ref) {
          call SplitMultiallelicVcf.SplitMultiallelicVcf as split_multiallelic_ref{
              input:
                  vcf_gz       = extract_chrom_variants_ref.chr_vcf,
                  vcf_idx      = extract_chrom_variants_ref.chr_vcf_idx,
                  output_prefix = "~{ref_prefix}.~{chromosomes[index]}.split",
                  docker_image = sv_pipeline_base_docker,
                  runtime_attr_override = runtime_attr_split_multiallelic
          }
        }

        File query_vcf_for_benchmark = select_first([
          split_multiallelic_query.output_vcf_gz,
          extract_chrom_variants_query.chr_vcf
        ])
        File query_vcf_idx_for_benchmark = select_first([
          split_multiallelic_query.output_vcf_idx,
          extract_chrom_variants_query.chr_vcf_idx
        ])
        File ref_vcf_for_benchmark = select_first([
          split_multiallelic_ref.output_vcf_gz,
          extract_chrom_variants_ref.chr_vcf
        ])
        File ref_vcf_idx_for_benchmark = select_first([
          split_multiallelic_ref.output_vcf_idx,
          extract_chrom_variants_ref.chr_vcf_idx
        ])

        if (defined(query_filter_vcf) && defined(ref_filter_vcf)){

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query_filter_a{
              input:
                  input_vcf = query_filter_vcf,
                  input_vcf_index = query_filter_vcf_idx,
                  chromosome = chromosomes[index],
                  output_name = "~{query_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                  docker_image = sv_pipeline_base_docker  
            }

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref_filter_a{
            input:
                input_vcf = ref_filter_vcf,
                input_vcf_index = ref_filter_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }  
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_query_ref{
              input:
                query_vcf = query_vcf_for_benchmark,
                query_vcf_idx = query_vcf_idx_for_benchmark,

                ref_vcf = ref_vcf_for_benchmark,
                ref_vcf_idx = ref_vcf_idx_for_benchmark,

                query_filter_vcf = extract_chrom_variants_query_filter_a.chr_vcf, 
                query_filter_vcf_idx = extract_chrom_variants_query_filter_a.chr_vcf_idx,

                ref_filter_vcf = extract_chrom_variants_ref_filter_a.chr_vcf, 
                ref_filter_vcf_idx = extract_chrom_variants_ref_filter_a.chr_vcf_idx, 

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }
        
        }

        if (defined(query_filter_vcf) && !defined(ref_filter_vcf)){
          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query_filter_b{
              input:
                  input_vcf = query_filter_vcf,
                  input_vcf_index = query_filter_vcf_idx,
                  chromosome = chromosomes[index],
                  output_name = "~{query_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                  docker_image = sv_pipeline_base_docker  
            }

          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_query{
              input:
                query_vcf = query_vcf_for_benchmark,
                query_vcf_idx = query_vcf_idx_for_benchmark,

                ref_vcf = ref_vcf_for_benchmark,
                ref_vcf_idx = ref_vcf_idx_for_benchmark,

                query_filter_vcf = extract_chrom_variants_query_filter_b.chr_vcf, 
                query_filter_vcf_idx = extract_chrom_variants_query_filter_b.chr_vcf_idx,

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }

        }

        if (!defined(query_filter_vcf) && defined(ref_filter_vcf)){

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref_filter_b{
            input:
                input_vcf = ref_filter_vcf,
                input_vcf_index = ref_filter_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }  
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_ref{
              input:
                query_vcf = query_vcf_for_benchmark,
                query_vcf_idx = query_vcf_idx_for_benchmark,

                ref_vcf = ref_vcf_for_benchmark,
                ref_vcf_idx = ref_vcf_idx_for_benchmark,

                ref_filter_vcf = extract_chrom_variants_ref_filter_b.chr_vcf, 
                ref_filter_vcf_idx = extract_chrom_variants_ref_filter_b.chr_vcf_idx, 

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }
        
        }

        if (!defined(query_filter_vcf) && !defined(ref_filter_vcf)){
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig{
              input:
                query_vcf = query_vcf_for_benchmark,
                query_vcf_idx = query_vcf_idx_for_benchmark,

                ref_vcf = ref_vcf_for_benchmark,
                ref_vcf_idx = ref_vcf_idx_for_benchmark,

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }

        }


      Array[File] tp_query_list = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.tp_query,benchmark_individual_vcf_per_contig_filter_query.tp_query,benchmark_individual_vcf_per_contig_filter_ref.tp_query,benchmark_individual_vcf_per_contig.tp_query])
      Array[File] tp_ref_list = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.tp_ref,benchmark_individual_vcf_per_contig_filter_query.tp_ref,benchmark_individual_vcf_per_contig_filter_ref.tp_ref,benchmark_individual_vcf_per_contig.tp_ref])
      Array[File] fp_query_list = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.fp_query,benchmark_individual_vcf_per_contig_filter_query.fp_query,benchmark_individual_vcf_per_contig_filter_ref.fp_query,benchmark_individual_vcf_per_contig.fp_query])
      Array[File] fp_ref_list = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.fp_ref,benchmark_individual_vcf_per_contig_filter_query.fp_ref,benchmark_individual_vcf_per_contig_filter_ref.fp_ref,benchmark_individual_vcf_per_contig.fp_ref])
    }

    scatter (i in range(length(sample_ids))){
        call ExtractFileByIndex.ExtractFileByIndex as extract_tp_query{
          input:
            index = i,
            nested_files = tp_query_list
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_query{
          input:
            vcfs = extract_tp_query.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.tp_query",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_tp_ref{
          input:
            index = i,
            nested_files = tp_ref_list
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_ref{
          input:
            vcfs = extract_tp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.tp_ref",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_fp_query{
          input:
            index = i,
            nested_files = fp_query_list
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_query{
          input:
            vcfs = extract_fp_query.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.fp_query",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_fp_ref{
          input:
            index = i,
            nested_files = fp_ref_list
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_ref{
          input:
            vcfs = extract_fp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.fp_ref",
            sv_base_mini_docker = sv_base_mini_docker
        }

        # --- merge TP + FP calls per sample into one tagged VCF, for query and ref ---
        call MergeTagVcf as merge_tag_query {
          input:
            tp_vcf = combine_vcfs_tp_query.concat_vcf,
            fp_vcf = combine_vcfs_fp_query.concat_vcf,
            sample_id = sample_ids[i],
            source_label = "~{query_prefix}.vs.~{ref_prefix}.query",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call MergeTagVcf as merge_tag_ref {
          input:
            tp_vcf = combine_vcfs_tp_ref.concat_vcf,
            fp_vcf = combine_vcfs_fp_ref.concat_vcf,
            sample_id = sample_ids[i],
            source_label = "~{query_prefix}.vs.~{ref_prefix}.ref",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call LongReadGenotypeTasks.CalcuCompStat as calcu_comp_stat{
          input:
            tp_query = combine_vcfs_tp_query.concat_vcf,
            tp_ref = combine_vcfs_tp_ref.concat_vcf,
            fp_query = combine_vcfs_fp_query.concat_vcf,
            fp_ref = combine_vcfs_fp_ref.concat_vcf,
            docker_image = sv_pipeline_base_docker
        }

        call LongReadGenotypeTasks.PlotCompResults as plot_comp_results{
          input:
            comp_stat = calcu_comp_stat.comp_stat,
            docker_image = sv_pipeline_base_docker
        }
    }

    # --- cohort-level FPR by chromosome, split by genomic context, for query and ref ---
    call CalcAndPlotFPRByContext as fpr_query {
      input:
        merged_vcfs = merge_tag_query.merged_vcf,
        label = "query",
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call CalcAndPlotFPRByContext as fpr_ref {
      input:
        merged_vcfs = merge_tag_ref.merged_vcf,
        label = "ref",
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }

  output {
    Array[File] merged_query_vcf = merge_tag_query.merged_vcf
    Array[File] merged_query_vcf_idx = merge_tag_query.merged_vcf_idx
    Array[File] merged_ref_vcf = merge_tag_ref.merged_vcf
    Array[File] merged_ref_vcf_idx = merge_tag_ref.merged_vcf_idx

    File fpr_table_query = fpr_query.fpr_table
    File fpr_plot_query = fpr_query.fpr_plot
    File fpr_table_ref = fpr_ref.fpr_table
    File fpr_plot_ref = fpr_ref.fpr_plot

    Array[File] benchmark_stat = calcu_comp_stat.comp_stat
    Array[File] benchmark_figure = plot_comp_results.figure
  }

}

task MergeTagVcf {
    input {
        File tp_vcf
        File fp_vcf
        String sample_id
        String source_label
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = "~{sample_id}.~{source_label}.merged"

    command <<<
        set -euo pipefail

        # tag TP records
        zcat ~{tp_vcf} | awk 'BEGIN{OFS="\t"}
            /^##/ {print; next}
            /^#CHROM/ {
                print "##INFO=<ID=benchmark,Number=1,Type=String,Description=\"Benchmark classification: TP or FP\">"
                print
                next
            }
            {
                if ($8 == ".") { $8 = "benchmark=TP" } else { $8 = $8";benchmark=TP" }
                print
            }' | bgzip -c > tp.tagged.vcf.gz

        # tag FP records (same header addition so headers match for concat)
        zcat ~{fp_vcf} | awk 'BEGIN{OFS="\t"}
            /^##/ {print; next}
            /^#CHROM/ {
                print "##INFO=<ID=benchmark,Number=1,Type=String,Description=\"Benchmark classification: TP or FP\">"
                print
                next
            }
            {
                if ($8 == ".") { $8 = "benchmark=FP" } else { $8 = $8";benchmark=FP" }
                print
            }' | bgzip -c > fp.tagged.vcf.gz

        tabix -p vcf tp.tagged.vcf.gz
        tabix -p vcf fp.tagged.vcf.gz

        bcftools concat -a tp.tagged.vcf.gz fp.tagged.vcf.gz \
          | bcftools sort -Oz -o ~{output_prefix}.vcf.gz
        tabix -p vcf ~{output_prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
    }
}

task CalcAndPlotFPRByContext {
    input {
        Array[File] merged_vcfs
        String label
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # extract CHROM, GC context annotation, and benchmark tag from every sample's merged VCF
        while read -r vcf; do
            bcftools query -f '%CHROM\t%INFO/GC\t%INFO/benchmark\n' "$vcf" >> combined.~{label}.tsv
        done < ~{write_lines(merged_vcfs)}

        cat <<'EOF' > run_fpr.R
        args <- commandArgs(trailingOnly = TRUE)
        in_file  <- args[1]
        out_table <- args[2]
        out_plot  <- args[3]
        label     <- args[4]

        df <- read.table(in_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                          col.names = c("chrom", "GC", "benchmark"))

        df$context_group <- NA
        df$context_group[df$GC %in% c("US", "RM")] <- "GC_US_RM"
        df$context_group[df$GC %in% c("SD", "SR")] <- "GC_SD_SR"
        df <- df[!is.na(df$context_group), ]

        chrom_order <- paste0("chr", c(1:22, "X", "Y"))
        df$chrom <- factor(df$chrom, levels = chrom_order[chrom_order %in% unique(df$chrom)])

        res <- do.call(rbind, lapply(split(df, list(df$chrom, df$context_group)), function(sub) {
          if (nrow(sub) == 0) return(NULL)
          data.frame(
            chrom = as.character(unique(sub$chrom)),
            context_group = unique(sub$context_group),
            n_total = nrow(sub),
            n_FP = sum(sub$benchmark == "FP"),
            FPR = sum(sub$benchmark == "FP") / nrow(sub)
          )
        }))
        res$chrom <- factor(res$chrom, levels = chrom_order[chrom_order %in% unique(res$chrom)])
        res <- res[order(res$context_group, res$chrom), ]
        write.csv(res, out_table, row.names = FALSE)

        chroms <- levels(res$chrom)
        groupA <- res[res$context_group == "GC_US_RM", ]
        groupB <- res[res$context_group == "GC_SD_SR", ]
        groupA <- groupA[match(chroms, groupA$chrom), ]
        groupB <- groupB[match(chroms, groupB$chrom), ]

        pdf(out_plot, width = 9, height = 5)
        ymax <- max(c(groupA$FPR, groupB$FPR), na.rm = TRUE) * 1.2
        plot(seq_along(chroms), groupA$FPR, type = "b", col = "#2E5C8A", pch = 16, lwd = 2,
             xaxt = "n", ylim = c(0, ymax),
             xlab = "Chromosome", ylab = "FPR = FP / (TP+FP)",
             main = paste0("False positive rate by chromosome and genomic context (", label, ")"))
        axis(1, at = seq_along(chroms), labels = chroms, las = 2, cex.axis = 0.7)
        lines(seq_along(chroms), groupB$FPR, type = "b", col = "#B23A48", pch = 17, lwd = 2)
        legend("topright",
               legend = c("GC=US/RM (unique / repeat-masked)", "GC=SD/SR (segdup / simple repeat)"),
               col = c("#2E5C8A", "#B23A48"), pch = c(16, 17), lwd = 2, bty = "n", cex = 0.8)
        dev.off()
EOF

        Rscript run_fpr.R combined.~{label}.tsv fpr_table.~{label}.csv fpr_plot.~{label}.pdf ~{label}
    >>>

    output {
        File fpr_table = "fpr_table.~{label}.csv"
        File fpr_plot = "fpr_plot.~{label}.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 7.5,
        disk_gb: 30,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_pipeline_base_docker
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}
