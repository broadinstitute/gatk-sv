version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectLargeCNVSupportForCPX {
    input {
        File cpx_ctx_bed
        String prefix

        File sample_batch_pe_map
        Array[String] batch_name_list
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_extract_cpx_lg_cnv_by_batch
        RuntimeAttr? runtime_attr_seek_depth_supp_for_cpx
        RuntimeAttr? runtime_attr_concat_bed_Step1
        RuntimeAttr? runtime_attr_concat_bed_Step2

    }

    call GenerateCnvSegmentFromCpx {
        input:
            bed = cpx_ctx_bed, #input file including cpx calls in bed format; 
            sample_batch_pe_map = sample_batch_pe_map,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_generate_cnv_segments_from_cpx
    }

    scatter (i in range(length(batch_name_list))){
        call ExtractCpxLgCnvByBatch {
            input:
                bed_gz = GenerateCnvSegmentFromCpx.cpx_lg_cnv,
                batch_name = batch_name_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_extract_cpx_lg_cnv_by_batch
        }

        call SeekDepthSuppForCpx as seek_depth_supp_for_cpx_del {
            input:
                cpx_lg_cnv = ExtractCpxLgCnvByBatch.lg_cnv_del,
                raw_depth_bed = Depth_DEL_beds[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_seek_depth_supp_for_cpx
        }

        call SeekDepthSuppForCpx as seek_depth_supp_for_cpx_dup {
            input:
                cpx_lg_cnv = ExtractCpxLgCnvByBatch.lg_cnv_dup,
                raw_depth_bed = Depth_DUP_beds[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_seek_depth_supp_for_cpx
        }

        call MiniTasks.ConcatBeds as concat_beds_svtype {
            input:
                shard_bed_files = [seek_depth_supp_for_cpx_del.cpx_cnv_depth_supp, seek_depth_supp_for_cpx_dup.cpx_cnv_depth_supp],
                prefix = "~{batch_name_list[i]}.lg_cnv.depth_supp",
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_concat_bed_Step1
        }
    }


    call MiniTasks.ConcatBeds as concat_beds_across_batches {
        input:
            shard_bed_files = concat_beds_svtype.merged_bed_file,
            prefix = "~{prefix}.lg_cnv.depth_supp",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_bed_Step2
        }

    output {
        File lg_cnv_depth_supp = concat_beds_across_batches.merged_bed_file
    }
}

task SeekDepthSuppForCpx {
    input {
        File cpx_lg_cnv
        File raw_depth_bed
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

    String prefix = basename(cpx_lg_cnv, '.bed.gz')

    output {
           File cpx_cnv_depth_supp = "~{prefix}.depth_supp.bed.gz"
    }
    command <<<

        set -e pipefail

        zcat ~{cpx_lg_cnv} | awk '{print $6}' | sort | uniq > sample_list.tsv

        echo -e '#chr\tpos\tend\tSVTYPE\tSVID\tsample\tbatch\tdepth_cov' > ~{prefix}.depth_supp.bed

        while read sample_name; do
            zcat ~{cpx_lg_cnv}    | awk -v sample="$sample_name" '{if ($6==sample) print}'  > query.bed
            zcat ~{raw_depth_bed} | awk -v sample="$sample_name" '{if ($5==sample) print}'  > ref.bed
            bedtools coverage -a query.bed -b ref.bed \
                | awk '{print $1,$2,$3,$4,$5,$6,$7,$NF}' \
                | sed -e 's/ /\t/g' \
                >> ~{prefix}.depth_supp.bed
        done < sample_list.tsv

        bgzip ~{prefix}.depth_supp.bed

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


task ExtractCpxLgCnvByBatch {
    input {
        File bed_gz
        String batch_name
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
           File lg_cnv_del = "~{batch_name}.lg_cnv.DEL.bed.gz"
           File lg_cnv_dup = "~{batch_name}.lg_cnv.DUP.bed.gz"
    }
    command <<<

        set -e pipefail

        zcat ~{bed_gz} | awk '{if ($4=="DEL" && $7=="~{batch_name}") print}' > ~{batch_name}.lg_cnv.DEL.bed
        zcat ~{bed_gz} | awk '{if ($4=="DUP" && $7=="~{batch_name}") print}' > ~{batch_name}.lg_cnv.DUP.bed
        bgzip ~{batch_name}.lg_cnv.DEL.bed
        bgzip ~{batch_name}.lg_cnv.DUP.bed

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

task GenerateCnvSegmentFromCpx {
    input {
        File bed
        File sample_batch_pe_map
        String sv_pipeline_docker
        Int min_depth_size = 5000
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(bed, ".bed.gz")

    command <<<
        set -euo pipefail

        python <<CODE

        import os
        import sys
        def split_cpx_interval(info):
            #eg of info: DUP_chrY:3125606-3125667
            out = [info.split('_')[0], info.split('_')[1].split(':')[0]] +  [int(i) for i in info.split('_')[1].split(':')[1].split('-')]
            return out[1:4]+[out[0]]

        def get_sample_batch_map(sample_batch_pe_map):
            sample_to_batch = {}
            with open(sample_batch_pe_map, 'r') as inp:
                for line in inp:
                    sample, batch, pe_file = line.strip().split("\t")
                    sample_to_batch[sample] = batch
            return sample_to_batch

        def readin_cpx_cnv(bed_input,sample_to_batch,min_depth_size):
            CPX_CNV = []
            f_bed = os.popen(r'''zcat %s'''%(bed_input))
            for line in f_bed:
                pin=line.strip().split('\t')
                if pin[0][0]=='#':
                    pos_CPX_TYPE = pin.index('CPX_TYPE')
                    pos_CPX_INTERVALS = pin.index('CPX_INTERVALS')
                    pos_SOURCE = pin.index('SOURCE')
                    pos_SAMPLES = pin.index('samples')
                else:
                    interval = []
                    if 'DEL_' in pin[pos_CPX_INTERVALS] or 'DUP_' in pin[pos_CPX_INTERVALS]:
                        interval += [split_cpx_interval(i)+[pin[3]] for i in pin[pos_CPX_INTERVALS].split(",") if "DEL_" in i or "DUP_" in i]
                    if 'DEL_' in pin[pos_SOURCE] or 'DUP_' in pin[pos_SOURCE]:
                        interval += [split_cpx_interval(i)+[pin[3]] for i in pin[pos_SOURCE].split(",") if "DEL_" in i or "DUP_" in i]
                    if len(interval)>0:
                        for i in interval:
                            if i[2]-i[1]>=min_depth_size:
                                sample_names = pin[pos_SAMPLES].split(',')
                                if i[3]=="DEL":
                                    for j in sample_names:
                                        CPX_CNV.append(i+[j, sample_to_batch[j]])
                                if i[3]=="DUP":
                                    for j in sample_names:
                                        CPX_CNV.append(i+[j, sample_to_batch[j]])
            f_bed.close()
            return CPX_CNV

        def write_cpx_cnv(CPX_CNV, fileout):
            fo=open(fileout, 'w')
            for info in CPX_CNV:
                print('\t'.join([str(i) for i in info]), file=fo)
            fo.close()

        bed_input = "~{bed}"
        fileout = "~{prefix}.lg_CNV.bed"
        min_depth_size = ~{min_depth_size}
        sample_to_batch = get_sample_batch_map("~{sample_batch_pe_map}")
        CPX_CNV = readin_cpx_cnv(bed_input, sample_to_batch, min_depth_size)
        write_cpx_cnv(CPX_CNV, fileout)
        CODE

        bgzip "~{prefix}.lg_CNV.bed"
    >>>

    output {
        File cpx_lg_cnv = "~{prefix}.lg_CNV.bed.gz"
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

