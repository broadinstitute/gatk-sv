version 1.0

import "Structs.wdl"
import "GATKSVTools.wdl" as gatk

workflow GATKSVPreprocessSample {
  input {
    # Sample data (must provide at least one file)
    String sample_id
    File? manta_vcf
    File? melt_vcf
    File? wham_vcf
    File? delly_vcf
    File cnv_bed
    # Note this will not work on gCNV VCFs generated with GATK < v4.1.5.0 ! Use cnv_beds instead
    File? gcnv_segments_vcf

    # Filtering options
    Int min_svsize = 50
    Int gcnv_min_qs = 80

    # Defragmentation options
    Float defrag_padding_fraction = 0.25

    # Reference
    File ref_dict
    File ref_fasta_fai

    # Docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_base_docker

    # VM resource options
    RuntimeAttr? runtime_attr_standardize
    RuntimeAttr? runtime_attr_defrag
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  Array[File?] optional_sv_vcfs_ = [manta_vcf, melt_vcf, wham_vcf, delly_vcf]
  Array[String] optional_sv_algorithms_ = ["manta", "melt", "wham", "delly"]
  scatter (i in range(length(optional_sv_algorithms_))) {
    if (defined(optional_sv_vcfs_[i])) {
      File scattered_sv_vcfs_ = select_first([optional_sv_vcfs_[i]])
      String scattered_sv_algorithms_ = optional_sv_algorithms_[i]
    }
  }
  Array[File] sv_vcfs_ = select_all(scattered_sv_vcfs_)
  Array[String] sv_algorithms_ = select_all(scattered_sv_algorithms_)

  call StandardizeVcfs {
    input:
      vcfs = sv_vcfs_,
      algorithms = sv_algorithms_,
      sample = sample_id,
      gcnv_segments_vcf = gcnv_segments_vcf,
      cnv_bed = cnv_bed,
      output_basename = "~{sample_id}.preprocess_sample",
      gcnv_min_qs = gcnv_min_qs,
      min_svsize = min_svsize,
      ref_fasta_fai = ref_fasta_fai,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_standardize
  }

  call gatk.SVCluster as DefragmentCNVs {
    input:
      vcfs = [StandardizeVcfs.out_cnv],
      vcf_indexes = [StandardizeVcfs.out_cnv_index],
      output_name="~{sample_id}.defragmented",
      ref_dict=ref_dict,
      vid_prefix="~{sample_id}_cnv_",
      omit_members=true,
      fast_mode=true,
      depth_sample_overlap=0,
      defrag_padding_fraction=defrag_padding_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_defrag
  }

  call gatk.SVCluster as MergeSVs {
    input:
      vcfs = [StandardizeVcfs.out_pesr, DefragmentCNVs.out],
      vcf_indexes = [StandardizeVcfs.out_pesr_index, DefragmentCNVs.out_index],
      output_name="~{sample_id}.std",
      ref_dict=ref_dict,
      vid_prefix="~{sample_id}_",
      omit_members=true,
      fast_mode=true,
      depth_overlap_fraction=1,
      mixed_overlap_fraction=1,
      pesr_overlap_fraction=1,
      depth_breakend_window=0,
      mixed_breakend_window=0,
      pesr_breakend_window=0,
      depth_sample_overlap=0,
      mixed_sample_overlap=0,
      pesr_sample_overlap=0,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_pesr
  }

  output {
    File out = MergeSVs.out
    File out_index = MergeSVs.out_index
  }
}

task StandardizeVcfs {
  input {
    Array[File] vcfs
    Array[String] algorithms
    File? gcnv_segments_vcf
    File? cnv_bed
    String sample
    String output_basename
    Int gcnv_min_qs
    Int min_svsize
    File ref_fasta_fai
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_vcfs = length(vcfs)

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 3
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euxo pipefail

    ############################################################
    # Filter and standardize gCNV segments VCF
    ############################################################

    if ~{defined(gcnv_segments_vcf)}; then
      # Force GT type to String to avoid a bcftools bug
      tabix ~{gcnv_segments_vcf}
      zcat ~{gcnv_segments_vcf} \
        | sed 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' \
        | bgzip \
        > gcnv_reheadered.vcf.gz
      tabix gcnv_reheadered.vcf.gz

      # With older gCNV versions, piping vcf directly into bcftools gives this error:
      # [W::vcf_parse] Contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)
      # Note use of --no-version to avoid header timestamp, which breaks call caching
      bcftools view \
        --no-version \
        -s ~{sample} \
        -i 'ALT!="." && QS>~{gcnv_min_qs}' \
        -O z \
        -o gcnv_filtered.vcf.gz \
        gcnv_reheadered.vcf.gz
      tabix gcnv_filtered.vcf.gz

      # Standardize by adding required INFO fields
      # Note this will not work on VCFs generated with GATK < v4.1.5.0 !
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/END\t%ID\n' gcnv_filtered.vcf.gz \
        | awk -F "\t" -v OFS="\t" '{
          if ($4=="<DEL>")  {
            svtype="DEL"; strands="+-";
          } else if ($4=="<DUP>") {
            svtype="DUP"; strands="-+";
          } else {
            svtype="."; strands=".";
          }
          print $1,$2,$3,$4,$5,$6,$5-$2+1,"depth",svtype,strands
        }' \
        | bgzip \
        > ann.tab.gz
      tabix -s1 -b2 -e2 ann.tab.gz

      echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">' > header_lines.txt
      echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> header_lines.txt
      echo '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness [++,+-,-+,--]">' >> header_lines.txt
      echo '##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">' >> header_lines.txt

      bcftools annotate \
        --no-version \
        -a ann.tab.gz \
        -h header_lines.txt \
        -c CHROM,POS,REF,ALT,INFO/END,ID,INFO/SVLEN,INFO/ALGORITHMS,INFO/SVTYPE,INFO/STRANDS \
        -O z \
        -o gcnv.vcf.gz \
        gcnv_filtered.vcf.gz
      tabix gcnv.vcf.gz
      echo "gcnv.vcf.gz" > cnv_vcfs.list
    fi

    ############################################################
    # Standardize SV caller VCFs
    ############################################################

    vcfs=(~{sep=" " vcfs})
    algorithms=(~{sep=" " algorithms})
    for (( i=0; i<~{num_vcfs}; i++ )); do
      vcf=${vcfs[$i]}
      algorithm=${algorithms[$i]}
      tabix $vcf
      svtk standardize \
        --sample-names ~{sample} \
        --prefix ${algorithm}_~{sample} \
        --contigs ~{ref_fasta_fai} \
        --min-size ~{min_svsize} \
        $vcf unsorted.vcf ${algorithm}
      bcftools sort unsorted.vcf -O z -o sorted.${algorithm}.vcf.gz
      tabix sorted.${algorithm}.vcf.gz
      echo "sorted.${algorithm}.vcf.gz" >> vcfs.list
    done

    ############################################################
    # Convert bed files into a standardized VCF
    ############################################################

    if ~{defined(cnv_bed)}; then
      echo "~{sample}" > samples.list
      # filter, concat, and header a combined bed
      zcat ~{cnv_bed} \
        | awk -F "\t" -v OFS="\t" '{if ($5=="~{sample}") print}' \
        | grep -v ^# \
        | sort -k1,1V -k2,2n \
        > cnvs.bed

      # note svtk generates an index automatically
      svtk rdtest2vcf --contigs ~{ref_fasta_fai} cnvs.bed samples.list cnvs.vcf.gz
      echo "cnvs.vcf.gz" >> cnv_vcfs.list
    fi

    ############################################################
    # Combine and sanitize VCFs
    ############################################################

    bcftools concat --no-version -a --file-list vcfs.list \
      | bcftools annotate --no-version -x ^FORMAT/GT \
      | grep -v "^##GATKCommandLine=" \
      | awk -F '\t' -v OFS='\t' '{ if ($0~/^#/) {print; next;} if ($8~"SVTYPE=BND") {$5="<BND>"} $4="N"; $6="."; $7="."; print}' \
      | sed -E 's/<INS:.+>/<INS>/g' \
      | bgzip \
      > ~{output_basename}.pesr.vcf.gz
    tabix ~{output_basename}.pesr.vcf.gz

    bcftools concat --no-version -a --file-list cnv_vcfs.list \
      | bcftools annotate --no-version -x ^FORMAT/GT,INFO/AC,INFO/AN,INFO/CHR2,INFO/END2 \
      | grep -v "^##GATKCommandLine=" \
      | grep -v "^##source=" \
      | bgzip \
      > ~{output_basename}.cnv.vcf.gz
    tabix ~{output_basename}.cnv.vcf.gz
  >>>

  output {
    File out_pesr = "~{output_basename}.pesr.vcf.gz"
    File out_pesr_index = "~{output_basename}.pesr.vcf.gz.tbi"
    File out_cnv = "~{output_basename}.cnv.vcf.gz"
    File out_cnv_index = "~{output_basename}.cnv.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
