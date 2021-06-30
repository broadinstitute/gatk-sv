version 1.0

import "Structs.wdl"

workflow DropRedundantCNVs {
  input {
    File vcf
    String contig
    String sv_pipeline_docker
  }

  call DropRedundantCNVs_1 {
    input:
      vcf=vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_2 {
    input:
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_3 {
    input:
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_4 {
    input:
      intervals_preclustered_subset_bed=DropRedundantCNVs_2.intervals_preclustered_subset_bed,
      step2_intervals_preclustered_subset_txt=DropRedundantCNVs_3.step2_intervals_preclustered_subset_txt,
      samples_list=DropRedundantCNVs_1.samples_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_5 {
    input:
      vcf=vcf,
      vids_to_remove_list_1=DropRedundantCNVs_4.vids_to_remove_list_1,
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      step2_variants_to_resolve_list=DropRedundantCNVs_4.step2_variants_to_resolve_list,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_6 {
    input:
      unsorted_vcf=DropRedundantCNVs_5.unsorted_vcf,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File cleaned_vcf_shard = DropRedundantCNVs_6.cleaned_vcf_shard
  }
}

task DropRedundantCNVs_1 {
  input {
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    bcftools query --list-samples ~{vcf} > samples.list

    ###PREP FILES
    #Convert full VCF to BED intervals
    #Ignore CPX events with UNRESOLVED filter status
    svtk vcf2bed --split-cpx --info SVTYPE \
      <(bcftools view -e 'INFO/SVTYPE == "CPX" && FILTER == "UNRESOLVED"' ~{vcf}) out.bed
    grep -e '^#\|DEL\|DUP\|CNV\|CPX' out.bed \
      | awk -v OFS="\t" '{ if ($5=="CN0") print $1, $2, $3, $4, "DEL", $5"\n"$1, $2, $3, $4, "DUP", $5; \
        else if ($5=="DEL" || $5=="DUP") print $1, $2, $3, $4, $6, $5 }' \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      | bgzip -c \
      > intervals.preclustered.bed.gz
  >>>

  output {
    File intervals_preclustered_bed = "intervals.preclustered.bed.gz"
    File samples_list = "samples.list"
  }
}

task DropRedundantCNVs_2 {
  input {
    File intervals_preclustered_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###REMOVE CNVS REDUNDANT WITH COMPLEX EVENTS
    #Subset to only variants that share some overlap (at least 10% recip) with at least one CPX variant
    bedtools intersect -wa -r -f 0.1 \
      -a ~{intervals_preclustered_bed} \
      -b <( zcat ~{intervals_preclustered_bed} | fgrep "CPX" ) \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      | uniq \
      | bgzip -c \
      > intervals.preclustered.subset.bed.gz
  >>>

  output {
    File intervals_preclustered_subset_bed = "intervals.preclustered.subset.bed.gz"
  }
}


task DropRedundantCNVs_3 {
  input {
    File intervals_preclustered_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 20.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###FIND REMAINING REDUNDANT CNVS WITH STRONG (80%) OVERLAP IN SAMPLES AND SIZE
    #Find CNV intervals that have 80% reciprocal overlap
    bedtools intersect -wa -wb -r -f 0.8 \
      -a ~{intervals_preclustered_bed} \
      -b ~{intervals_preclustered_bed} \
      | awk -v FS="\t" '{ if ($4!=$10 && $6==$12) print $0 }' \
      | awk -v OFS="\t" '$4 ~ /DEL|DUP/ { print $0 }' \
      | awk -v OFS="\t" '$10 ~ /DEL|DUP/ { print $0 }' \
      | cut -f4,5,10,11 \
      | sort \
      | uniq \
      | gzip \
      > step2.intervals.preclustered.subset.txt.gz
  >>>

  output {
    File step2_intervals_preclustered_subset_txt = "step2.intervals.preclustered.subset.txt.gz"
  }
}


task DropRedundantCNVs_4 {
  input {
    File step2_intervals_preclustered_subset_txt
    File intervals_preclustered_subset_bed
    File samples_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(step2_intervals_preclustered_subset_txt, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 60,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE | bgzip > out.vcf.gz
    import sys
    import gzip
    from collections import namedtuple, defaultdict

    import pysam
    import numpy as np
    from scipy import sparse
    from scipy.sparse import csgraph

    BedCall = namedtuple('BedCall', 'chrom start end name samples svtype'.split())

    def reciprocal_overlap(a, b, frac):
      if a.chrom != b.chrom:
        return False
      if a.start >= b.end or b.start >= a.end:
        return False
      ov = min(a.end, b.end) - max(a.start, b.start)
      return (ov / float(max(a.end - a.start, b.end - b.start))) >= frac


    def sample_overlap(samples_a, samples_b, denom, frac):
      if len(samples_a) == 0 or len(samples_b) == 0:
        return True
      ov = len(samples_a.intersection(samples_b))
      return (ov / float(denom)) >= frac


    def read_intervals(path, samples_dict):
      intervals = []
      with gzip.open(path, "rb") as f:
        for lineb in f:
          tokens = lineb.decode('utf-8').strip().split('\t')
          sample_indexes = set([samples_dict[s] for s in tokens[4].split(',')])
          intervals.append(BedCall(tokens[0], int(tokens[1]), int(tokens[2]), tokens[3], sample_indexes, tokens[5]))
        return intervals

    # Save memory using sample id indexing
    with open("~{samples_list}") as f:
      samples_list = [line.strip() for line in f]
    num_samples = len(samples_list)
    samples_dict = {samples_list[i]: i for i in range(num_samples)}

    intervals = read_intervals("~{intervals_preclustered_subset_bed}", samples_dict)
    num_intervals = len(intervals)

    # 50% RO and sample overlap in subsetted intervals
    # Generate sparse graph for clustering
    RO_FRAC = 0.5
    G = sparse.eye(len(intervals), dtype=np.uint8, format='lil')
    for i in range(num_intervals):
      ro_indexes = [j for j in range(i) if reciprocal_overlap(intervals[i], intervals[j], RO_FRAC)]
      for j in ro_indexes:
        G[i, j] = 1

    # Compute clusters
    n_comp, cluster_labels = csgraph.connected_components(G, connection='weak', directed=False)
    clusters = defaultdict(list)
    for i in range(len(cluster_labels)):
      clusters[cluster_labels[i]].append(i)

    # Find CNVs in clusters containing at least one CPX
    SAMPLE_FRAC = 0.5
    vids_to_remove = set([])
    for cluster in clusters.values():
      cnvs = [i for i in cluster if "DEL" in intervals[i].name or "DUP" in intervals[i].name]
      cpx = [i for i in cluster if "CPX" in intervals[i].name]
      for i in cnvs:
        for j in cpx:
          if sample_overlap(intervals[i].samples, intervals[j].samples, len(intervals[i].samples), SAMPLE_FRAC):
            vids_to_remove.add(intervals[i].name + "\n")
            break

    with open("VIDs_to_remove.list", 'w') as f:
      f.writelines(sorted(list(vids_to_remove)))

    # Find clusters of CNVs only, using 80% overlap parameters
    with gzip.open("~{step2_intervals_preclustered_subset_txt}") as f:
      intervals2 = []
      for line in f:
        tokens = line.decode('utf-8').strip().split('\t')
        samples_a = set([samples_dict[s] for s in tokens[1].split(',')])
        samples_b = set([samples_dict[s] for s in tokens[3].split(',')])
        intervals2.append((tokens[0], samples_a, tokens[2], samples_b))

    num_intervals2 = len(intervals2)
    vids_to_resolve_list = []
    SAMPLE_FRAC2 = 0.8
    for interval in intervals2:
      samples_a = interval[1]
      samples_b = interval[3]
      union = samples_a.union(samples_b)
      if sample_overlap(samples_a, samples_b, len(union), SAMPLE_FRAC2):
        vids_to_resolve_list.append("{}\n".format(",".join(sorted([interval[0], interval[2]]))))

    vids_to_resolve_list = sorted(list(set(vids_to_resolve_list)))

    with open("step2.variants_to_resolve.list", 'w') as f:
      f.writelines(vids_to_resolve_list)

    CODE
  >>>

  output {
    File step2_variants_to_resolve_list = "step2.variants_to_resolve.list"
    File vids_to_remove_list_1 = "VIDs_to_remove.list"
  }
}


task DropRedundantCNVs_5 {
  input {
    File vcf
    File vids_to_remove_list_1
    File intervals_preclustered_bed
    File step2_variants_to_resolve_list
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, intervals_preclustered_bed, intervals_preclustered_bed, step2_variants_to_resolve_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 30,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    python3 <<CODE | bgzip > drop_redundant_cnvs_5.~{contig}.vcf.gz
    import sys
    import pysam
    import gzip

    sys.stderr.write("Reading step2...\n")
    with open("~{step2_variants_to_resolve_list}") as f:
        vids_sets_to_resolve = [set(line.strip().split(',')) for line in f.readlines()]
        vids_list = sorted(list(set([x for y in vids_sets_to_resolve for x in y])))

    sys.stderr.write("Reading vids to remove...\n")
    with open("~{vids_to_remove_list_1}") as f:
        vids_to_remove = set([line.strip() for line in f.readlines()])

    sys.stderr.write("Reading preclustered intervals...\n")
    with gzip.open("~{intervals_preclustered_bed}") as f:
        intervals = {}
        for lineb in f:
            tokens = lineb.decode('utf-8').strip().split('\t')
            vid = tokens[3]
            intervals[vid] = tokens

    sys.stderr.write("Finding partners...\n")
    partners = {}
    all_partners = set([])
    for vid in vids_list:
        # get all other variants from clusters containing this variant
        partners[vid] = set([p for vset in vids_sets_to_resolve if vid in vset for p in vset])
        all_partners.update(partners[vid])

    vids_to_remove.update(all_partners)
    #with open("vids_to_remove_2.list", 'w') as f:
    #    f.writelines(sorted([v+"\n" for v in vids_to_remove]))

    sys.stderr.write("Scanning vcf...\n")
    vcf = pysam.VariantFile("~{vcf}")
    records = {r.id: r for r in vcf if r.id in all_partners}
    vcf.close()

    def count_gts(record):
        result = [0, 0, 0]
        num_samples = len(record.samples)
        for g in [record.samples[i]['GT'] for i in range(num_samples)]:
            if g == (0, 0):
                result[1] += 1
            elif g == (None, None):
                result[2] += 1
            else:
                result[0] += 1
        return result

    def get_best_score_vid(scores):
        return sorted(scores.items(), key=lambda x: x[1])[-1][0]

    sys.stderr.write("Generating records...\n")
    records_to_add = []
    processed_vids = set([])
    for vid in vids_list:
        if vid in processed_vids:
            continue
        vid_partners = partners[vid]
        processed_vids.update(vid_partners)
        partner_intervals = [intervals[p] for p in vid_partners]
        most_samples_vid = sorted(partner_intervals, key=lambda x : len(x[4].split(',')))[-1][3]
        x = sorted(partner_intervals, key=lambda x : len(x[4].split(',')))
        best_genotype_vid = None
        best_non_ref = -1
        best_ref = -1
        scores = {p: count_gts(records[p]) for p in vid_partners}
        scores_non_ref = {p: scores[p][0] for p in vid_partners if scores[p][0] > 0}
        scores_ref = {p: scores[p][1] for p in vid_partners if scores[p][1] > 0}
        scores_no_call = {p: scores[p][2] for p in vid_partners if scores[p][2] > 0}
        if len(scores_non_ref) > 0:
            best_genotype_vid = get_best_score_vid(scores_non_ref)
        elif len(scores_ref) > 0:
            best_genotype_vid = get_best_score_vid(scores_ref)
        else:
            best_genotype_vid = get_best_score_vid(scores_no_call)
        sys.stderr.write(most_samples_vid + "\n")
        s1 = str(records[most_samples_vid]).split('\t')[0:9]
        s2 = str(records[best_genotype_vid]).split('\t', 9)
        records_to_add.append("\t".join(s1) + "\t" + s2[9])

    sys.stderr.write("Writing vcf...\n")
    vcf = pysam.VariantFile("~{vcf}")
    sys.stdout.write(str(vcf.header))
    for record in vcf:
      if record.id not in vids_to_remove:
        sys.stdout.write(str(record))
    vcf.close()

    for record in records_to_add:
      sys.stdout.write(record)

    CODE

  >>>

  output {
    File unsorted_vcf = "drop_redundant_cnvs_5.~{contig}.vcf.gz"
  }
}


task DropRedundantCNVs_6 {
  input {
    File unsorted_vcf
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = contig + ".shard.no_CNV_redundancies.vcf.gz"

  Float input_size = size(unsorted_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 20.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###CLEAN UP FINAL OUTPUT
    zcat ~{unsorted_vcf} \
      | vcf-sort \
      | bgzip \
      > ~{outfile_name}
  >>>

  output {
    File cleaned_vcf_shard = outfile_name
  }
}

