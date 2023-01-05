version 1.0

import "Structs.wdl"

# Get lists of PCRPLUS and PCRMINUS samples present in input VCF
task GetSampleLists {
  input{
    String sv_base_mini_docker
    File vcf
    File vcf_idx
    File? pcrplus_samples_list
    String prefix
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10 + size(vcf, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    bcftools query -l ~{vcf} > all_samples.list
    if ~{defined(pcrplus_samples_list)}; then
      awk -v OFS="\t" 'ARGIND==1{inFileA[$1]; next} {if($1 in inFileA){print $1,"PCRPLUS"}else{print $1,"PCRMINUS"}}' ~{pcrplus_samples_list} all_samples.list \
        > ~{prefix}.PCR_status_assignments.txt
    else
      awk -v OFS="\t" '{ print $1, "PCRMINUS" }' all_samples.list \
        > ~{prefix}.PCR_status_assignments.txt
    fi
  >>>

  output {
    File sample_PCR_labels = "~{prefix}.PCR_status_assignments.txt"
  }

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

# Restrict famfile to probands from complete trios with all members present in the VCF
task SubsetFamfile {
  input {
    File famfile
    File sample_PCR_labels
    String sv_base_mini_docker
  }

  runtime {
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 20 HDD"
    bootDiskSizeGb: 10
    docker: sv_base_mini_docker
    preemptible: 3
    maxRetries: 1
  }

  command <<<
    set -eu -o pipefail

    # Get list of samples from families not present in VCF
    cut -f1 ~{sample_PCR_labels} | sort | uniq > all_samples.list
    cut -f2-4 ~{famfile} | sed 's/\t/\n/g' | sort | uniq > famfile_samples.list
    fgrep -wvf all_samples.list famfile_samples.list > missing_samples.list

    # Restrict famfile to proband entries from trios without missing samples
    fgrep -wvf missing_samples.list ~{famfile} \
    | awk -v FS="\t" -v OFS="\t" \
      '{ if ($2!="0" && $3!="0" && $4!="0") print $0 }' \
    > cleaned_trios_in_vcf.fam
  >>>

  output {
    File subsetted_famfile = "cleaned_trios_in_vcf.fam"
  }
}

# Concatenate minGQ tarballs trained on different chromosomes
task ConcatTarball {
  input {
    Array[File] tarballs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2, 
    mem_gb: 3, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -e

    while read tarball; do
      tar xzf $tarball
    done < ~{write_lines(tarballs)}
    find . -name trio_variant_info.txt.gz > files.txt
    mkdir -p out
    NUM_SHARDS=$(cat files.txt | sed 's/\/cacheCopy//g' | sed 's/\/attempt-[2-9]//g' | tr '/' '\t' | cut -f7 | sort | uniq | wc -l)

    for (( i=0; i<$NUM_SHARDS; i++ )); do
      fgrep "/shard-$i/" files.txt > shard.txt
      mkdir -p "out/shard-$i/"
      OUT="out/shard-$i/trio_variant_info.txt.gz"
      sed 1q shard.txt | xargs -n1 gunzip -c | bgzip > $OUT
      cat shard.txt | tail -n+2 | xargs -n1 gunzip -c | { grep -v ^# || true; } | bgzip >> $OUT
    done

    echo "test1"
    tar -czf "out.tar.gz" out/
  >>>

  output {
    File tarball = "out.tar.gz"
  }

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

# Downsample minGQ tarballs to a random set of N trios
task SubsetTrioTarball {
  input {
    File tarball_in
    Int max_trios
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2, 
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu -o pipefail

    mkdir tar_contents/
    tar -xzvf ~{tarball_in} --directory tar_contents/

    mkdir ~{prefix}/
    find tar_contents/ -type d -name "shard-*" \
    | shuf --random-source=<( yes 2022 ) \
    | head -n ~{max_trios} \
    | xargs -I {} mv {} ~{prefix}/
    tar -czvf ~{prefix}.tar.gz ~{prefix}
  >>>

  output {
    File subsetted_tarball = "~{prefix}.tar.gz"
  }

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

# Split a VCF into two parts, corresponding to PCR+ and PCR-
task SplitPcrVcf {
  input{
    File vcf
    String prefix
    File? pcrplus_samples_list
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    if [ ! -z "~{pcrplus_samples_list}" ] && [ $( cat "~{pcrplus_samples_list}" | wc -l ) -gt 0 ]; then
      #Get index of PCR+ samples
      PCRPLUS_idxs=$( zcat ~{vcf} | sed -n '1,500p' | fgrep "#" | fgrep -v "##" \
                      | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' \
                      | fgrep -wf ~{pcrplus_samples_list} | cut -f1 | paste -s -d, )
      #Get PCR+ VCF
      zcat ~{vcf} \
      | cut -f1-9,"$PCRPLUS_idxs" \
      | bgzip -c \
      > "~{prefix}.PCRPLUS.vcf.gz"
      tabix -f -p vcf "~{prefix}.PCRPLUS.vcf.gz"
      #Get PCR- VCF
      zcat ~{vcf} \
      | cut --complement -f"$PCRPLUS_idxs" \
      | bgzip -c \
      > "~{prefix}.PCRMINUS.vcf.gz"
      tabix -f -p vcf "~{prefix}.PCRMINUS.vcf.gz"
    else
      cp ~{vcf} ~{prefix}.PCRMINUS.vcf.gz
      tabix -f -p vcf "~{prefix}.PCRMINUS.vcf.gz"
      touch ~{prefix}.PCRPLUS.vcf.gz
      touch ~{prefix}.PCRPLUS.vcf.gz.tbi
    fi
  >>>

  output {
    File PCRPLUS_vcf = "~{prefix}.PCRPLUS.vcf.gz"
    File PCRPLUS_vcf_idx = "~{prefix}.PCRPLUS.vcf.gz.tbi"
    File PCRMINUS_vcf = "~{prefix}.PCRMINUS.vcf.gz"
    File PCRMINUS_vcf_idx = "~{prefix}.PCRMINUS.vcf.gz.tbi"
  }

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


# Get a simple table with ID/AC/AN/AF per variant, prior to minGQ
task GetAfTables {
  input{
    File vcf
    File vcf_idx
    String prefix
    String sv_pipeline_docker
    File? pcrplus_samples_list
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10 + size(vcf, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    #Run vcf2bed
    svtk vcf2bed --info ALL --no-samples ~{vcf} "~{prefix}.vcf2bed.bed"
    #Cut to necessary columns
    idxs=$( sed -n '1p' "~{prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print $1, NR }' \
            | grep -e 'name\|SVLEN\|SVTYPE\|_AC\|_AN\|_CN_NONREF_COUNT\|_CN_NUMBER' \
            | fgrep -v "OTH" \
            | cut -f2 \
            | paste -s -d\, || true )
    cut -f"$idxs" "~{prefix}.vcf2bed.bed" \
      | sed 's/^name/\#VID/g' \
      | gzip -c \
      > "~{prefix}.frequencies.preclean.txt.gz"
    if [ ! -z "~{pcrplus_samples_list}" ]; then
      echo -e "dummy\tPCRMINUS\ndummy2\tPCRPLUS" > dummy.tsv
    else
      echo -e "dummy\tPCRMINUS" > dummy.tsv
    fi
    #Clean frequencies
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/clean_frequencies_table.R \
      "~{prefix}.frequencies.preclean.txt.gz" \
      dummy.tsv \
      "~{prefix}.frequencies.txt"
    for PCR in $( cut -f2 dummy.tsv | sort | uniq ); do
      AC_idx=$( zcat "~{prefix}.frequencies.txt.gz" | sed -n '1p' | sed 's/\t/\n/g' | awk -v PCR="$PCR" '{ if ($1==PCR"_AC") print NR }' )
      AN_idx=$( zcat "~{prefix}.frequencies.txt.gz" | sed -n '1p' | sed 's/\t/\n/g' | awk -v PCR="$PCR" '{ if ($1==PCR"_AN") print NR }' )
      zcat "~{prefix}.frequencies.txt.gz" \
        | sed '1d' \
        | awk -v FS="\t" -v OFS="\t" -v AC="$AC_idx" -v AN="$AN_idx" \
          '{ print $1, $(AC), $(AN) }' \
      > ~{prefix}."$PCR".AF_preMinGQ.txt
    done
    if [ ! -z ~{prefix}.PCRPLUS.AF_preMinGQ.txt ]; then
      touch ~{prefix}.PCRPLUS.AF_preMinGQ.txt
    fi
  >>>

  output {
    File PCRPLUS_AF_table = "~{prefix}.PCRPLUS.AF_preMinGQ.txt"
    File PCRMINUS_AF_table = "~{prefix}.PCRMINUS.AF_preMinGQ.txt"
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

# Shard a trio famfile to keep only trios that are all represented in the vcf header
task SplitFamfile {
  input{
    File? vcf
    File? vcf_idx
    Int? max_count_famfile_shards
    File famfile
    String prefix
    String sv_pipeline_docker
    Int fams_per_shard
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    #Get list of sample IDs & column numbers from VCF header
    tabix -H ~{vcf} | fgrep -v "##" | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt

    set -euo pipefail
    python3 <<CODE > "~{prefix}.cleaned_trios.fam"
    import random
    import sys
    random.seed(928483)
    SAMPLES_PATH="vcf_header_columns.txt"
    PED_PATH="~{famfile}"
    SHARD_SIZE=int(~{fams_per_shard})
    MAX_SHARDS=int(~{max_count_famfile_shards})

    with open(SAMPLES_PATH) as f:
      samples = set([line.strip().split()[0] for line in f])
    with open(PED_PATH) as f:
      ped = []
      for line in f:
        family = line.strip().split('\t') 
        if family[1] in samples and family[2] in samples and family[3] in samples:
          ped.append(family)

    if len(ped) > SHARD_SIZE * MAX_SHARDS:
      random.shuffle(ped)
      ped = ped[:SHARD_SIZE * MAX_SHARDS]

    for family in ped:
      sys.stdout.write("\t".join(family) + "\n")
    CODE

    split -l ~{fams_per_shard} --numeric-suffixes=00001 -a 5 ~{prefix}.cleaned_trios.fam famfile_shard_
  >>>

  output {
    File cleaned_trios_famfile = "~{prefix}.cleaned_trios.fam"
    Array[File] famfile_shards = glob("famfile_shard_*")
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

# Collect a single table of all relevant variants for a single family
task CollectTrioSVdat {
  input{
    Array[File] vcf_shard_idxs
    File famfile
    String filter_metric = "GQ"
    Array[String] gather_trio_geno_options = []
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2, 
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    # Make list of all samples in famfile
    awk -v OFS="\n" '{ print $2, $3, $4 }' ~{famfile} > samples.list

    for wrapper in 1; do
      #Write header
      echo -e "#famID\tVID\tSVLEN\tAF\tSVTYPE\tFILTER\tpro_EV\tpro_AC\tfa_AC\tmo_AC\tpro_~{filter_metric}\tfa_~{filter_metric}\tmo_~{filter_metric}"

      #Iterate over list of VCF shards
      while read vcf_idx; do

        # Localize VCF from remote bucket while subsetting to samples in trios
        vcf="$( echo $vcf_idx | sed 's/\.tbi$//g' | sed 's/\/cromwell_root/gs\:\//g' )"
        subset_vcf="$( basename $vcf | sed 's/\.vcf\.gz/\.subsetted\.vcf\.gz/g' )"
        mv $vcf_idx ./ #VCF index must be in execution directory for remote streaming to work
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools view \
          -O z \
          --output-file "$subset_vcf" \
          --samples-file samples.list \
          --min-ac 1 \
          --no-update \
          $vcf

        #Get list of sample IDs & column numbers from VCF header
        zfgrep "#" $subset_vcf | fgrep -v "##" | head -n1000 | sed 's/\t/\n/g' \
        | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt

        #Iterate over families & subset VCF
        while read famID pro fa mo prosex pheno; do

          #Subset vcf to only non-ref or multiallelic sites in the family
          bcftools view \
            -O z \
            --output-file $famID.vcf.gz \
            --samples "$pro","$fa","$mo" \
            --min-ac 1 \
            --exclude 'FILTER == "MULTIALLELIC"' \
            --no-update \
            "$subset_vcf"

          #Get list of CNVs in proband that are ≥5kb have ≥50% coverage in either parent
          svtk vcf2bed -i SVTYPE --no-header $famID.vcf.gz stdout \
          | awk -v OFS="\t" '{ if ($NF ~ /DEL|DUP|CNV/) print $1, $2, $3, $4, $NF, $6 }' \
          > $famID.CNVs.bed
          fgrep -w $pro $famID.CNVs.bed \
          | awk -v OFS="\t" '{ if ($3-$2>=5000 && $5!="CNV") print $1, $2, $3, $4, $5 }' \
          > $pro.CNVs.gt5kb.bed
          fgrep -w $fa $famID.CNVs.bed > $fa.CNVs.bed
          fgrep -w $mo $famID.CNVs.bed > $mo.CNVs.bed

          # Find deletions where AC needs to be updated due to parental coverage
          awk -v OFS="\t" '{ if ($NF=="DEL") print $0, "1" }' $pro.CNVs.gt5kb.bed \
          | bedtools coverage -a - \
            -b <( awk '{ if ($5 ~ /DEL|CNV/) print $0 }' $fa.CNVs.bed ) \
          | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
          | bedtools coverage -a - \
            -b <( awk '{ if ($5 ~ /DEL|CNV/) print $0 }' $mo.CNVs.bed ) \
          | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
          > $famID.RD_genotype_update.txt

          # Find duplications where AC needs to be updated due to parental coverage
          awk -v OFS="\t" '{ if ($NF=="DUP") print $0, "1" }' $pro.CNVs.gt5kb.bed \
          | bedtools coverage -a - \
            -b <( awk '{ if ($5 ~ /DUP|CNV/) print $0 }' $fa.CNVs.bed ) \
          | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
          | bedtools coverage -a - \
            -b <( awk '{ if ($5 ~ /DUP|CNV/) print $0 }' $mo.CNVs.bed ) \
          | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
          >> $famID.RD_genotype_update.txt

          # Get variant stats
          /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/gather_trio_genos.py \
            --ac-adj $famID.RD_genotype_update.txt \
            --metric ~{filter_metric} \
            ~{sep=" " gather_trio_geno_options} \
            --no-header \
            $famID.vcf.gz stdout "$pro" "$fa" "$mo" \
          | awk -v famID="$famID" -v OFS="\t" '{ print famID, $0 }'

        done < ~{famfile}
      done < ~{write_lines(vcf_shard_idxs)}
    done | gzip -c > "trio_variant_info.txt.gz"
  >>>

  output {
    File trio_SVdata = "trio_variant_info.txt.gz"
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


# Gather all trio SV data into a single tarball (helps with Cromwell file localization)
task GatherTrioData {
  input{
    Array[File] files
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 3 * ceil(size(files, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    tar -czvf ~{prefix}.tar.gz -T ~{write_lines(files)}
  >>>

  output {
    File tarball = "~{prefix}.tar.gz"
  }

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


# Enumerate all minGQ conditions to test
task EnumerateConditions {
  input{
    String prefix
    Int condition_shards
    String optimize_minSizes
    String optimize_maxSizes
    String optimize_minFreqs
    String optimize_maxFreqs
    String optimize_includeSVTYPEs
    String optimize_includeFILTERs
    String optimize_excludeFILTERs
    String optimize_includeEV
    String optimize_excludeEV
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

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_tranches_table.R \
      --min.sizes "~{optimize_minSizes}" \
      --max.sizes "~{optimize_maxSizes}" \
      --min.freqs "~{optimize_minFreqs}" \
      --max.freqs "~{optimize_maxFreqs}" \
      --svtype.include "~{optimize_includeSVTYPEs}" \
      --filter.include "~{optimize_includeFILTERs}" \
      --filter.exclude "~{optimize_excludeFILTERs}" \
      --ev.include "~{optimize_includeEV}" \
      --ev.exclude "~{optimize_excludeEV}" \
      "~{prefix}.minGQ_conditions.txt"
    fgrep -v "#" "~{prefix}.minGQ_conditions.txt" \
      > "~{prefix}.minGQ_conditions.noHeader.txt"
    /opt/sv-pipeline/04_variant_resolution/scripts/evenSplitter.R \
      -S ~{condition_shards} \
      --shuffle \
      "~{prefix}.minGQ_conditions.noHeader.txt" \
      "~{prefix}.minGQ_conditions.noHeader.shard"
  >>>

  output {
    File conditions_table = "~{prefix}.minGQ_conditions.txt"
    File conditions_table_noHeader = "~{prefix}.minGQ_conditions.noHeader.txt"
    Array[File] conditions_table_noHeader_shards = glob("~{prefix}.minGQ_conditions.noHeader.shard*")
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


# Merge ROC optimal cutoffs or stats
task CombineRocOptResults {
  input{
    Array[String] shards
    String outfile
    String gcloud_sdk_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Note: shards is intentionally localized as Array[String] and not Array[File] to avoid long (multi-hour) call cache wait times
    # Its a little clunky, but its faster and works fine
    cat ~{write_lines(shards)} \
    | xargs -I {} gsutil cat {} \
    | fgrep -v "#" \
    | sort -Vk1,1 \
    | uniq \
    > ~{outfile}
  >>>

  output {
    File merged_file = "~{outfile}"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gcloud_sdk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


# Build final minGQ filtering tree
task BuildFilterTree {
  input{
    File conditions_table
    File condition_optimizations
    File condition_distrib_stats
    String prefix
    String optimize_metric = "GQ"
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_lookup_table.R \
      --optimize-metric "~{optimize_metric}" \
      "~{conditions_table}" \
      "~{condition_distrib_stats}" \
      "~{condition_optimizations}" \
      "~{prefix}.min~{optimize_metric}.ordered_tree_hierarchy.txt" \
      "~{prefix}.min~{optimize_metric}.filter_lookup_table.txt"
  >>>


  output {
    File ordered_tree_hierarchy = "~{prefix}.min~{optimize_metric}.ordered_tree_hierarchy.txt"
    File filter_lookup_table = "~{prefix}.min~{optimize_metric}.filter_lookup_table.txt"
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


# Apply minMetric filter to VCF
task FilterGTs {
  input{
    File vcf
    File? minMetric_lookup_table
    String prefix
    String PCR_status
    Float maxNCR
    Array[String] filter_GT_options = []
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override    
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu -o pipefail

    cmd="/opt/sv-pipeline/scripts/downstream_analysis_and_filtering/filter_GTs_by_metric.py "
    cmd="$cmd ~{sep=' ' filter_GT_options} "
    cmd="$cmd --max-ncr ~{maxNCR} --simplify-INS-SVTYPEs --cleanAFinfo "
    cmd="$cmd --prefix \"~{PCR_status}\" ~{vcf} ~{minMetric_lookup_table} stdout"

    echo -e "Filtering GTs with the following command:\n$cmd"
    
    eval $cmd \
    | fgrep -v "##INFO=<ID=AN," \
    | fgrep -v "##INFO=<ID=AC," \
    | fgrep -v "##INFO=<ID=AF," \
    | fgrep -v "##INFO=<ID=N_BI_GENOS," \
    | fgrep -v "##INFO=<ID=N_HOMREF," \
    | fgrep -v "##INFO=<ID=N_HET," \
    | fgrep -v "##INFO=<ID=N_HOMALT," \
    | fgrep -v "##INFO=<ID=FREQ_HOMREF," \
    | fgrep -v "##INFO=<ID=FREQ_HET," \
    | fgrep -v "##INFO=<ID=FREQ_HOMALT," \
    | fgrep -v "##INFO=<ID=CN_NUMBER," \
    | fgrep -v "##INFO=<ID=CN_COUNT," \
    | fgrep -v "##INFO=<ID=CN_FREQ," \
    | fgrep -v "##INFO=<ID=CN_NONREF_COUNT," \
    | fgrep -v "##INFO=<ID=CN_NONREF_FREQ," \
    | bgzip -c \
    > "~{prefix}.filtered.vcf.gz"
    tabix -p vcf -f "~{prefix}.filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "~{prefix}.filtered.vcf.gz"
    File filtered_vcf_idx = "~{prefix}.filtered.vcf.gz.tbi"
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


# Merge PCRPLUS and PCRMINUS VCFs for a single chromosome
task MergePcrVCFs {
  input{
    File? PCRPLUS_vcf
    File PCRMINUS_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    if [ ! -z "~{PCRPLUS_vcf}" ];then
      #Sanitize FILTER columns
      zcat "~{PCRPLUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRPLUS_filters.txt
      zcat "~{PCRMINUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/merge_filter_columns.py \
        PCRPLUS_filters.txt \
        PCRMINUS_filters.txt \
        merged_filters.txt
      #Write new VCF header
      zgrep -e '^##' ~{PCRPLUS_vcf} > "~{prefix}.minGQ_filtered.vcf"
      zgrep -e '^##' ~{PCRMINUS_vcf} | fgrep "NOCALL_RATE" >> "~{prefix}.minGQ_filtered.vcf"
      #Column-wise merger
      paste \
        <( zgrep -ve '^##' "~{PCRPLUS_vcf}" | cut -f1-6 ) \
        <( cat <( echo -e "FILTER" ) merged_filters.txt ) \
        <( zgrep -ve '^##' "~{PCRPLUS_vcf}" | cut -f8- ) \
        <( zgrep -ve '^##' "~{PCRMINUS_vcf}"  | cut -f10- ) \
      >> "~{prefix}.minGQ_filtered.vcf"
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "~{prefix}.minGQ_filtered.vcf" \
        "~{prefix}.minGQ_filtered.no_blanks.vcf"
      #Bgzip & tabix
      bgzip -f "~{prefix}.minGQ_filtered.no_blanks.vcf"
      tabix -p vcf -f "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
    else
      #Sanitize FILTER columns
      zcat "~{PCRMINUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
      #Write new VCF header
      zcat "~{PCRMINUS_vcf}" | sed -n '1,1000p' | grep -e '^##' > "~{prefix}.minGQ_filtered.vcf"
      #Column-wise merger
      paste \
        <( zcat "~{PCRMINUS_vcf}" | grep -ve '^##' | cut -f1-6 ) \
        <( cat <( echo -e "FILTER" ) PCRMINUS_filters.txt ) \
        <( zcat "~{PCRMINUS_vcf}" | grep -ve '^##' | cut -f8- ) \
      >> "~{prefix}.minGQ_filtered.vcf"
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "~{prefix}.minGQ_filtered.vcf" \
        "~{prefix}.minGQ_filtered.no_blanks.vcf"
      #Bgzip & tabix
      vcf-sort "~{prefix}.minGQ_filtered.no_blanks.vcf" | bgzip > "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
      tabix -p vcf -f "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
    fi
  >>>

  output {
    File merged_vcf = "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
    File merged_vcf_idx = "~{prefix}.minGQ_filtered.no_blanks.vcf.gz.tbi"
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
