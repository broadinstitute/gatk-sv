version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RunVisualizePlots.wdl" as VisualizePlots

workflow ManualReview {
  input {
    File vcf
    File plot_stratifications
    String prefix
    
    # Optional sharding
    Int? variants_per_shard
    
    # Docker images
    String sv_pipeline_docker
    String igv_docker
    
    # Optional reference genes for LOEUF filtering
    File? reference_genes
    
    # Optional flags for RdTest
    String? flags
    
    # IGV-specific inputs (optional)
    Int? igv_max_window
    String? igv_buffer
    File? batch_bincov
    File? sample_batches
    File? batch_medianfile
    File? rd_outliers
    File? sample_crai_cram
    File? reference
    File? reference_index
    Boolean? file_localization
    Boolean? requester_pays
    String? sv_base_mini_docker
    String? variant_interpretation_docker
    
    # Depth plotting inputs (conditional based on stratifications)
    Array[File]? median_files
    Array[File]? rd_files
    File? ped_file
    Int? min_size_rd
    
    # Runtime attributes
    RuntimeAttr? runtime_attr_filter_vcf
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_igv
    RuntimeAttr? runtime_attr_aggregate_plots
    RuntimeAttr? runtime_attr_scatter_vcf
  }

  # Parse stratifications to determine what plot types are needed
  call ParseStratifications {
    input:
      plot_stratifications = plot_stratifications,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_filter_vcf
  }

  # Shard VCF if requested
  if (defined(variants_per_shard)) {
    call MiniTasks.ScatterVcf {
      input:
        vcf = vcf,
        prefix = prefix,
        records_per_shard = select_first([variants_per_shard]),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_scatter_vcf
    }
  }
  
  Array[File] vcf_shards = if defined(variants_per_shard) then select_first([ScatterVcf.shards]) else [vcf]

  # Filter VCF based on stratifications for each shard
  scatter (i in range(length(vcf_shards))) {
    call FilterVcfForPlots {
      input:
        vcf_file = vcf_shards[i],
        plot_stratifications = plot_stratifications,
        reference_genes = reference_genes,
        prefix = prefix + "_shard_" + i,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_filter_vcf
    }
  }

  # Flatten and aggregate all filtered VCFs by plot name
  call AggregateFilteredVcfs {
    input:
      filter_summaries = FilterVcfForPlots.summary_file,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker
  }

  # Create plots using appropriate method based on requirements
  if (ParseStratifications.needs_depth_plots || ParseStratifications.needs_igv_plots) {
    File ped_file_required = select_first([ped_file])
    
    # Use RunVisualizePlots for IGV stratifications (handles both RD and IGV)
    if (ParseStratifications.needs_igv_plots) {
      scatter (igv_plot_info in AggregateFilteredVcfs.igv_plot_info) {
        call CreateVarfile {
          input:
            vcf_file = igv_plot_info.vcf_file,
            plot_name = igv_plot_info.plot_name,
            sv_pipeline_docker = sv_pipeline_docker
        }
        
        # Convert rd_files/median_files to batch format if needed
        if (!defined(batch_bincov) && defined(rd_files)) {
          call CreateBatchFiles {
            input:
              rd_files = select_first([rd_files]),
              median_files = select_first([median_files]),
              prefix = igv_plot_info.plot_name,
              sv_pipeline_docker = sv_pipeline_docker
          }
        }
        
        call VisualizePlots.VisualizePlots as IgvVisualizePlots {
          input:
            varfile = CreateVarfile.varfile,
            pedfile = ped_file_required,
            prefix = igv_plot_info.plot_name,
            run_RD = true,
            run_IGV = true,
            run_evidence_plots = false,
            run_cram_plots = true,
            batch_bincov = select_first([batch_bincov, CreateBatchFiles.batch_bincov_file]),
            sample_batches = sample_batches,
            batch_medianfile = select_first([batch_medianfile, CreateBatchFiles.batch_medianfile_file]),
            rd_outliers = rd_outliers,
            sample_crai_cram = sample_crai_cram,
            buffer = select_first([igv_buffer, "1000"]),
            igv_max_window = igv_max_window,
            reference = reference,
            reference_index = reference_index,
            file_localization = file_localization,
            requester_pays = requester_pays,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_rdtest_docker = sv_pipeline_docker,
            igv_docker = igv_docker,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_rdtest = runtime_attr_igv
        }
      }
    }
    
    # Use simple CreateDepthPlots for depth-only stratifications
    if (ParseStratifications.needs_depth_plots && !ParseStratifications.needs_igv_plots) {
      Array[File] median_files_required = select_first([median_files])
      Array[File] rd_files_required = select_first([rd_files])
      
      scatter (depth_plot_info in AggregateFilteredVcfs.depth_plot_info) {
        call CreateDepthPlots {
          input:
            bed_file = depth_plot_info.bed_file,
            plot_name = depth_plot_info.plot_name,
            median_files = median_files_required,
            rd_files = rd_files_required,
            ped_file = ped_file_required,
            min_size = select_first([min_size_rd, 50]),
            flags = flags,
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_rdtest
        }
      }
    }
  }

  # Aggregate all plots into final tarball
  call AggregatePlots {
    input:
      depth_plot_tarballs = select_all(CreateDepthPlots.plots_tarball),
      igv_plot_tarballs = select_all(IgvVisualizePlots.output_plots),
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_aggregate_plots
  }

  output {
    File plots_tarball = AggregatePlots.plots_tarball
  }
}

task ParseStratifications {
  input {
    File plot_stratifications
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Check if depth plots are needed
    if grep -q "True" <(cut -f15 ~{plot_stratifications} | tail -n+2); then
      echo "true" > needs_depth_plots.txt
    else
      echo "false" > needs_depth_plots.txt
    fi
    
    # Check if IGV plots are needed
    if grep -q "True" <(cut -f16 ~{plot_stratifications} | tail -n+2); then
      echo "true" > needs_igv_plots.txt
    else
      echo "false" > needs_igv_plots.txt
    fi
  >>>

  output {
    Boolean needs_depth_plots = read_boolean("needs_depth_plots.txt")
    Boolean needs_igv_plots = read_boolean("needs_igv_plots.txt")
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

task FilterVcfForPlots {
  input {
    File vcf_file
    File plot_stratifications
    File? reference_genes
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf_file, plot_stratifications], "GB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: ceil(10 + input_size * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    python3 /opt/sv-pipeline/scripts/filter_vcf_for_plots.py \
      ~{vcf_file} \
      ~{plot_stratifications} \
      ~{prefix} \
      ~{if defined(reference_genes) then "--reference-genes " + reference_genes else ""}

  >>>

  output {
    File summary_file = "~{prefix}.summary.tsv"
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

task AggregateFilteredVcfs {
  input {
    Array[File] filter_summaries
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    python3 << 'EOF'
import csv
import json
import os

# Combine all filter summaries and group by plot name
depth_plots = {}
igv_plots = {}

for summary_file in [~{sep=", " filter_summaries}]:
    summary_file = summary_file.strip('"')
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            plot_name = row['plot_name']
            vcf_file = row['vcf_file']
            
            if row['depth_plot'] == 'True':
                if plot_name not in depth_plots:
                    depth_plots[plot_name] = []
                depth_plots[plot_name].append(vcf_file)
            elif row['igv_plot'] == 'True':
                if plot_name not in igv_plots:
                    igv_plots[plot_name] = []
                igv_plots[plot_name].append(vcf_file)

# Create merged VCFs and BED files for each plot name
depth_plot_info = []
igv_plot_info = []

for plot_name, vcf_files in depth_plots.items():
    # Merge VCFs and convert to BED for depth plotting
    bed_file = f"{plot_name}.bed"
    with open(bed_file, 'w') as bed_out:
        for vcf_file in vcf_files:
            if os.path.exists(vcf_file) and os.path.getsize(vcf_file) > 0:
                with open(vcf_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            fields = line.strip().split('\t')
                            if len(fields) >= 8:
                                chrom, pos, var_id, ref, alt, qual, filter_field, info_field = fields[:8]
                                # Extract SVTYPE and SVLEN from INFO
                                svtype = ""
                                svlen = 1000  # default
                                for info in info_field.split(';'):
                                    if info.startswith('SVTYPE='):
                                        svtype = info.split('=')[1]
                                    elif info.startswith('SVLEN='):
                                        try:
                                            svlen = abs(int(info.split('=')[1]))
                                        except:
                                            pass
                                # BED format: chr start end id samples svtype
                                bed_out.write(f"{chrom}\t{pos}\t{int(pos)+svlen}\t{var_id}\tsample1\t{svtype}\n")
    
    depth_plot_info.append({"plot_name": plot_name, "bed_file": bed_file})

for plot_name, vcf_files in igv_plots.items():
    # Merge VCFs for IGV plotting
    merged_vcf = f"{plot_name}.vcf"
    header_written = False
    with open(merged_vcf, 'w') as vcf_out:
        for vcf_file in vcf_files:
            if os.path.exists(vcf_file) and os.path.getsize(vcf_file) > 0:
                with open(vcf_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            if not header_written:
                                vcf_out.write(line)
                        else:
                            vcf_out.write(line)
                    header_written = True
    
    igv_plot_info.append({"plot_name": plot_name, "vcf_file": merged_vcf})

# Write output arrays
with open("depth_plot_info.json", 'w') as f:
    json.dump(depth_plot_info, f)

with open("igv_plot_info.json", 'w') as f:
    json.dump(igv_plot_info, f)
EOF
  >>>

  output {
    Array[Object] depth_plot_info = read_json("depth_plot_info.json")
    Array[Object] igv_plot_info = read_json("igv_plot_info.json")
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

task CreateDepthPlots {
  input {
    File bed_file
    String plot_name
    Array[File] median_files
    Array[File] rd_files
    File ped_file
    Int min_size
    String? flags
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(bed_file, "GB") + size(rd_files, "GB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 7.5,
    disk_gb: ceil(100 + input_size * 10),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Index rd_files if needed
    for rd_file in ~{sep=" " rd_files}; do
      if [[ ! -f "${rd_file}.tbi" ]]; then
        tabix -p bed "$rd_file"
      fi
    done
    
    mkdir -p plots/~{plot_name}
    
    # Check if BED file has content
    if [[ -s ~{bed_file} ]]; then
      # Create median file list
      echo ~{sep=" " median_files} | tr ' ' '\t' > median_file.txt
      
      # Run RdTest with native support for other SVTYPEs
      Rscript /opt/RdTest/RdTestV2.R \
        -b ~{bed_file} \
        -n ~{plot_name} \
        -x "$(dirname ~{rd_files[0]})" \
        -m median_file.txt \
        -f ~{ped_file} \
        -p TRUE \
        ~{flags}
      
      # Move plots to organized directory
      mv *.jpg plots/~{plot_name}/ 2>/dev/null || true
    else
      echo "No variants for plot: ~{plot_name}"
    fi
    
    # Create tarball
    tar -czf ~{plot_name}_depth_plots.tar.gz plots/
  >>>

  output {
    File plots_tarball = "~{plot_name}_depth_plots.tar.gz"
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

task AggregatePlots {
  input {
    Array[File] depth_plot_tarballs
    Array[File] igv_plot_tarballs
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    mkdir -p ~{prefix}_manual_review_plots
    
    # Aggregate depth plots
    depth_tarballs=(~{sep=" " depth_plot_tarballs})
    for tarball in "${depth_tarballs[@]}"; do
      if [[ -f "$tarball" ]]; then
        tar -xzf "$tarball"
        # Move extracted plots to final directory
        find plots/ -type d -mindepth 1 -exec mv {} ~{prefix}_manual_review_plots/ \; 2>/dev/null || true
        rm -rf plots/
      fi
    done
    
    # Aggregate IGV plots
    igv_tarballs=(~{sep=" " igv_plot_tarballs})
    for tarball in "${igv_tarballs[@]}"; do
      if [[ -f "$tarball" ]]; then
        tar -xzf "$tarball"
        # Move extracted plots to final directory, merging if directory exists
        find plots/ -type d -mindepth 1 | while read plotdir; do
          plot_name=$(basename "$plotdir")
          if [[ -d "~{prefix}_manual_review_plots/$plot_name" ]]; then
            mv "$plotdir"/* "~{prefix}_manual_review_plots/$plot_name/" 2>/dev/null || true
          else
            mv "$plotdir" "~{prefix}_manual_review_plots/"
          fi
        done
        rm -rf plots/
      fi
    done
    
    # Create final tarball
    tar -czf ~{prefix}_manual_review_plots.tar.gz ~{prefix}_manual_review_plots/
  >>>

  output {
    File plots_tarball = "~{prefix}_manual_review_plots.tar.gz"
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

task CreateVarfile {
  input {
    File vcf_file
    String plot_name
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Convert VCF to varfile format expected by RunVisualizePlots
    echo -e "chrom\tstart\tend\tname\tsvtype\tsample" > ~{plot_name}.varfile.bed
    
    if [[ -s ~{vcf_file} ]]; then
      grep -v "^#" ~{vcf_file} | while read line; do
        fields=($line)
        chrom=${fields[0]}
        pos=${fields[1]}
        id=${fields[2]}
        ref=${fields[3]}
        alt=${fields[4]}
        info=${fields[7]}
        
        # Extract SVTYPE and SVLEN from INFO
        svtype=$(echo "$info" | grep -oP 'SVTYPE=\K[^;]+' || echo "UNK")
        svlen=$(echo "$info" | grep -oP 'SVLEN=\K[^;]+' || echo "1000")
        
        # Calculate end position
        if [[ $svlen =~ ^-?[0-9]+$ ]]; then
          end=$((pos + ${svlen#-}))
        else
          end=$((pos + 1000))
        fi
        
        # For each sample with variant, create a row
        echo -e "${chrom}\t${pos}\t${end}\t${id}\t${svtype}\tsample1" >> ~{plot_name}.varfile.bed
      done
    fi
  >>>

  output {
    File varfile = "~{plot_name}.varfile.bed"
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

task CreateBatchFiles {
  input {
    Array[File] rd_files
    Array[File] median_files
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Create batch_bincov file (batch name -> rd file path)
    echo -e "batch\tpath" > ~{prefix}_batch_bincov.txt
    i=1
    for rd_file in ~{sep=" " rd_files}; do
      echo -e "batch${i}\t${rd_file}" >> ~{prefix}_batch_bincov.txt
      i=$((i+1))
    done
    
    # Create batch_medianfile by concatenating all median files horizontally
    # First get the header from the first file
    head -n 1 ~{median_files[0]} > ~{prefix}_batch_medianfile.txt
    
    # Then get the data row from the first file
    tail -n 1 ~{median_files[0]} > temp_data.txt
    
    # Append data from subsequent files (skip headers)
    for median_file in ~{sep=" " median_files}; do
      if [[ "$median_file" != "~{median_files[0]}" ]]; then
        # Skip header, get data, and append horizontally
        tail -n 1 "$median_file" | cut -f2- | paste temp_data.txt - > temp_combined.txt
        mv temp_combined.txt temp_data.txt
      fi
    done
    
    # Combine header and data
    cat temp_data.txt >> ~{prefix}_batch_medianfile.txt
    rm -f temp_data.txt temp_combined.txt
  >>>

  output {
    File batch_bincov_file = "~{prefix}_batch_bincov.txt"
    File batch_medianfile_file = "~{prefix}_batch_medianfile.txt"
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
