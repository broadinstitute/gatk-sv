version 1.0
    
import "Structs.wdl"
import "ReformatRawFiles.wdl" as raw
import "TasksMakeCohortVcf.wdl" as miniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo
import "Utils.wdl" as util

workflow DeNovoSV {

    input {
        File vcf_file
        File? vcf_index
        File ped_input
        File batch_raw_file
        File batch_depth_raw_file
        File batch_bincov_index

        Array[String] contigs
        String prefix

        File? family_ids_txt

        File python_config
        File genomic_disorder_input

        File exclude_regions
        File sample_batches

        Int records_per_shard

        String variant_interpretation_docker
        String sv_base_mini_docker
        String sv_pipeline_updates_docker
        String python_docker
        RuntimeAttr? runtime_attr_gd
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_raw_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_merge_bed
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_divide_by_chrom
        RuntimeAttr? runtime_attr_raw_reformat_bed
        RuntimeAttr? runtime_attr_merge_final_bed_files
        RuntimeAttr? runtime_attr_create_plots
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_clean_ped
        RuntimeAttr? runtime_attr_call_outliers
        RuntimeAttr? runtime_attr_get_batched_files
        RuntimeAttr? runtime_attr_subset_by_samples
        RuntimeAttr? runtime_attr_merge_gd
        RuntimeAttr? runtime_attr_batch_vcf
        RuntimeAttr? runtime_attr_merge
    }

    # family_ids_txt is a text file containg one family id per line.
    # If this file is given, subset all other input files to only include the necessary batches.
    if (defined(family_ids_txt)) {
        File family_ids_txt_ = select_first([family_ids_txt])
        call GetBatchedFiles {
            input:
                batch_raw_file = batch_raw_file,
                batch_depth_raw_file = batch_depth_raw_file,
                ped_input = ped_input,
                family_ids_txt = family_ids_txt_,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                python_docker=python_docker,
                runtime_attr_override = runtime_attr_get_batched_files
        }

        call util.SubsetVcfBySamplesList {
            input:
                vcf = vcf_file,
                vcf_idx = vcf_index,
                list_of_samples = GetBatchedFiles.samples,
                outfile_name = prefix,
                keep_af = true,
                remove_private_sites = false,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_by_samples
        }
    }

    # Makes a ped file of singletons, duos, and trios for input into the de novo script (only including families of interest)
    call CleanPed {
        input:
            ped_input = ped_input,
            vcf_input = select_first([SubsetVcfBySamplesList.vcf_subset, vcf_file]),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_clean_ped
    }

    # Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call raw.ReformatRawFiles as ReformatRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([GetBatchedFiles.batch_raw_files_list, batch_raw_file]),
            ped_input = CleanPed.cleaned_ped,
            depth = false,
            variant_interpretation_docker = variant_interpretation_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    # Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call raw.ReformatRawFiles as ReformatDepthRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([GetBatchedFiles.batch_depth_raw_files_list, batch_depth_raw_file]),
            ped_input = CleanPed.cleaned_ped,
            depth = true,
            variant_interpretation_docker = variant_interpretation_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }
    
    scatter (i in range(length(contigs))) {
        # Generates a list of genomic disorder regions in the vcf input as well as in the depth raw files
        call GetGenomicDisorders {
            input:
                genomic_disorder_input=genomic_disorder_input,
                ped = ped_input,
                vcf_file = select_first([SubsetVcfBySamplesList.vcf_subset, vcf_file]),
                depth_raw_file_proband = ReformatDepthRawFiles.reformatted_proband_raw_files[i],
                depth_raw_file_parents = ReformatDepthRawFiles.reformatted_parents_raw_files[i],
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_gd
        }

        # Splits vcf by chromosome
        call SubsetVcf {
            input:
                vcf_file = select_first([SubsetVcfBySamplesList.vcf_subset, vcf_file]),
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        # Shards vcf
        call miniTasks.ScatterVcf as ScatterVcf {
            input:
                vcf=SubsetVcf.vcf_output,
                prefix=prefix,
                records_per_shard=records_per_shard,
                sv_pipeline_docker=sv_pipeline_updates_docker,
                runtime_attr_override=runtime_override_shard_vcf
        }
    
        # Runs the de novo calling python script on each shard and outputs a per chromosome list of de novo SVs
        call runDeNovo.DeNovoSVsScatter as GetDeNovo {
            input:
                ped_input=CleanPed.cleaned_ped,
                vcf_files=ScatterVcf.shards,
                disorder_input=GetGenomicDisorders.gd_output_for_denovo,
                chromosome=contigs[i],
                raw_proband=ReformatRawFiles.reformatted_proband_raw_files[i],
                raw_parents=ReformatRawFiles.reformatted_parents_raw_files[i],
                raw_depth_proband=ReformatDepthRawFiles.reformatted_proband_raw_files[i],
                raw_depth_parents=ReformatDepthRawFiles.reformatted_parents_raw_files[i],
                exclude_regions = exclude_regions,
                sample_batches = sample_batches,
                batch_bincov_index = select_first([GetBatchedFiles.batch_bincov_index_subset, batch_bincov_index]),
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_denovo = runtime_attr_denovo,
                runtime_attr_vcf_to_bed = runtime_attr_vcf_to_bed
        }
    }

    # Merges the per chromosome final de novo SV outputs
    call MergeDenovoBedFiles {
        input:
            bed_files = GetDeNovo.per_chromosome_final_output_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_final_bed_files
    }

    # Outputs a final callset of de novo SVs as well as outlier de novo SV calls
    call CallOutliers {
        input:
            bed_file = MergeDenovoBedFiles.merged_denovo_output,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_call_outliers
    }

    # Generates plots for QC
    call CreatePlots {
        input:
            bed_file = CallOutliers.final_denovo_nonOutliers_output,
            ped_input = ped_input,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_create_plots
    }

    # Merges the genomic disorder region output from each chromosome to compile a list of genomic disorder regions
    call MergeGenomicDisorders {
        input:
            genomic_disorder_input=GetGenomicDisorders.gd_output_from_depth_raw_files,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_gd
    }

    output {
        File cleaned_ped = CleanPed.cleaned_ped
        File final_denovo_nonOutliers = CallOutliers.final_denovo_nonOutliers_output
        File final_denovo_outliers = CallOutliers.final_denovo_outliers_output
        File final_denovo_nonOutliers_plots = CreatePlots.output_plots
        Array [File] denovo_output_annotated = GetDeNovo.per_chromosome_annotation_output_file
        File gd_depth = MergeGenomicDisorders.gd_output_from_depth
        File gd_vcf = GetGenomicDisorders.gd_output_from_final_vcf[1]
    }
}

task SubsetVcf {

    input {
        File vcf_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File vcf_output = "~{chromosome}.vcf.gz"
        File no_header_vcf_output = "~{chromosome}.noheader.vcf.gz"
    }

    command <<<
        set -euxo pipefail

        bcftools index ~{vcf_file}
        bcftools view ~{vcf_file} --regions ~{chromosome} -O z -o  ~{chromosome}.vcf.gz
        bcftools view ~{chromosome}.vcf.gz | grep -v ^## | bgzip -c > ~{chromosome}.noheader.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task GetGenomicDisorders {

    input {
        File vcf_file
        File ped
        File depth_raw_file_proband
        File depth_raw_file_parents
        String chromosome
        File genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([vcf_file, ped, genomic_disorder_input, depth_raw_file_parents, depth_raw_file_proband]), "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File gd_output_from_final_vcf = "gd.variants.from.final.vcf.txt.gz"
        File gd_output_from_depth_raw_files = "~{chromosome}.gd.variants.in.depth.raw.files.txt.gz"
        File gd_output_for_denovo = "annotated.gd.variants.names.txt"
    }

    command <<<
        set -euxo pipefail

        sort -k1,1 -k2,2n ~{genomic_disorder_input} > sorted.genomic.txt
        bedtools intersect -wa -wb -f 0.3 -r -a ~{vcf_file} -b ~{genomic_disorder_input} | cut -f 3 |sort -u > annotated.gd.variants.names.txt
        
        echo "Done with first line"

        bedtools intersect -wa -wb -f 0.3 -r -a ~{vcf_file} -b ~{genomic_disorder_input} > gd.variants.from.final.vcf.txt
        bgzip gd.variants.from.final.vcf.txt

        echo "Done with GD from vcf"
        
        Rscript /src/denovo/create_per_sample_bed.R ~{genomic_disorder_input} unsorted.gd.per.sample.txt unsorted.gd.per.family.txt ~{ped} ~{chromosome}
        sort -k1,1 -k2,2n unsorted.gd.per.sample.txt > gd.per.sample.txt
        sort -k1,1 -k2,2n unsorted.gd.per.family.txt > gd.per.family.txt
        cat ~{depth_raw_file_parents} | gunzip | sort -k1,1 -k2,2n | bgzip -c > sorted.depth.parents.bed.gz
        cat ~{depth_raw_file_proband} | gunzip | sort -k1,1 -k2,2n | bgzip -c > sorted.depth.proband.bed.gz

        echo "Done with R script"

        bedtools intersect -wa -wb -f 0.3 -r -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.proband.txt
        bedtools intersect -wa -wb -f 0.3 -r -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.parents.txt
        
        echo "done with intersect in depth variants"

        bedtools coverage -wa -wb -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz | awk '{if ($NF>=0.30) print }' > ~{chromosome}.coverage.parents.txt
        bedtools coverage -wa -wb -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz | awk '{if ($NF>=0.30) print }' > ~{chromosome}.coverage.proband.txt

        echo "done with coverage in depth variants"

        cat ~{chromosome}.coverage.parents.txt ~{chromosome}.coverage.proband.txt > ~{chromosome}.coverage.txt

        echo "done with cat"

        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt
        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt

        echo "done with intersect no -r"

        cat ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt > ~{chromosome}.remove.txt

        echo "done with cat"

        bedtools intersect -v -wb -b ~{chromosome}.remove.txt -a ~{chromosome}.coverage.txt > ~{chromosome}.kept.coverage.txt

        echo "done with grep"

        cat ~{chromosome}.gd.variants.in.depth.raw.file.proband.txt ~{chromosome}.gd.variants.in.depth.raw.file.parents.txt ~{chromosome}.kept.coverage.txt > ~{chromosome}.gd.variants.in.depth.raw.files.txt
        bgzip ~{chromosome}.gd.variants.in.depth.raw.files.txt
        echo "done with cat"
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task MergeGenomicDisorders {

    input {
        Array[File] genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(genomic_disorder_input, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File gd_output_from_depth = "gd.raw.files.output.txt.gz"
    }

    command {
        set -euxo pipefail

        zcat ~{sep=" " genomic_disorder_input} > gd.raw.files.output.txt
        bgzip gd.raw.files.output.txt
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task MergeDenovoBedFiles {

    input {
        Array[File] bed_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_files, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + (bed_files_size) * 2.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File merged_denovo_output = "denovo.merged.bed.gz"
    }

    command {
        set -exuo pipefail

        zcat ~{bed_files[0]} | awk 'NR<=1' > denovo.merged.bed
        zcat ~{sep=" " bed_files} | grep -v ^chrom >> denovo.merged.bed
        bgzip denovo.merged.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task CallOutliers {

    input {
        File bed_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_file, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 16, # 3.75
        disk_gb: ceil(10 + (bed_files_size) * 2.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File final_denovo_nonOutliers_output = "final.denovo.merged.bed.gz"
        File final_denovo_outliers_output = "final.denovo.merged.outliers.bed.gz"
        File final_denovo_allSamples_output = "final.denovo.merged.allSamples.bed.gz"
    }

    command {
        set -exuo pipefail

        python /src/denovo/denovo_outliers.py --bed ~{bed_file}

        bgzip final.denovo.merged.bed
        bgzip final.denovo.merged.outliers.bed
        bgzip final.denovo.merged.allSamples.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task CreatePlots {

    input {
        File bed_file
        File ped_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([bed_file, ped_input]), "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 16, # 3.75
        disk_gb: ceil(10 + input_size * 1.2),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File output_plots = "output_plots.pdf"
    }

    command {
        set -exuo pipefail

        Rscript /src/denovo/denovo_sv_plots.R ~{bed_file} ~{ped_input} output_plots.pdf
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task CleanPed {

    input {
        File ped_input
        File vcf_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float ped_size = size(ped_input, "GB")
    Float vcf_size = size(vcf_input, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + vcf_size + ped_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File cleaned_ped = "subset_cleaned_ped.txt"
    }

    command {
        set -exuo pipefail

        # TODO: this script should get the name of the output file it generates as an input argument.
        # The output filename is currently hardcoded to be 'clean_ped.txt'.
        Rscript /src/denovo/clean_ped.R ~{ped_input}
        cut -f2 cleaned_ped.txt | awk 'NR > 1' > all_samples.txt
        bcftools query -l ~{vcf_input} > samples_to_include_in_ped.txt

        # Note: grep returns a non-success exit code (i.e., other than `0`) when it cannot find the
        # match in the following scripts. We do not expect it to find a match for every entry.
        # Hence, to avoid exit with unsuccessful code, we can either drop pipefail from above or use `|| true`.

        grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt || true
        grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt || true
        grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt || true
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task GetBatchedFiles {

    input {
        File batch_raw_file
        File batch_depth_raw_file
        File family_ids_txt
        File ped_input
        File sample_batches
        File batch_bincov_index
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([batch_raw_file, batch_depth_raw_file, ped_input, sample_batches, batch_bincov_index, family_ids_txt]), "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File batch_raw_files_list = "batch_raw_files_list.txt"
        File batch_depth_raw_files_list = "batch_depth_raw_files_list.txt"
        File batch_bincov_index_subset = "batch_bincov_index.txt"
        File samples = "samples.txt"
    }

    command {
        set -exuo pipefail

        if grep -q -w -f ~{family_ids_txt} ~{ped_input}; then
            grep -w -f ~{family_ids_txt} ~{ped_input} | cut -f2 | sort -u > samples.txt
        else
            echo "No matching family IDs from family_ids_txt found in ped_input file." >&2
            exit 1
        fi

        if grep -q -w -f samples.txt ~{sample_batches}; then
            grep -w -f samples.txt ~{sample_batches} | cut -f2 | sort -u > batches.txt
        else
            echo "No matching individual IDs found in the sample_batches file." >&2
            exit 1
        fi

        grep -w -f batches.txt ~{batch_bincov_index} > batch_bincov_index.txt
        grep -w -f batches.txt ~{batch_raw_file} > batch_raw_files_list.txt
        grep -w -f batches.txt ~{batch_depth_raw_file} > batch_depth_raw_files_list.txt
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: python_docker
    }
}