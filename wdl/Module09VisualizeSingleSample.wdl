version 1.0

import "Structs.wdl"
import "RdTestVisualization.wdl" as rdtest
import "IGVTrioPlotsAllSamples.wdl" as igv_trio
import "IGVGeneratePlotsAllSamples.wdl" as igv_individual

workflow Module09VisualizeSingleSample{
    input{
        File Fasta
        File Fasta_idx
        File Fasta_dict

        File varfile
        File pedfile
        String flags
        String prefix
        File batch_bincov
        File sample_batches

        Array[File] medianfile
        Array[String] sample_list
        Array[File] cram_list
        Array[File] crai_list

        String sv_base_mini_docker
        String sv_pipeline_rdtest_docker
        String igv_docker

        RuntimeAttr? runtime_attr_override
        RuntimeAttr? runtime_attr_concatinate
        RuntimeAttr? runtime_attr_rdtest
        }
    
    call rdtest.RdTestVisualization as RdTest{
        input:
            flags = flags,
            prefix = prefix,
            bed = varfile,
            pedfile = pedfile,
            medianfile = medianfile,
            batch_bincov=batch_bincov,
            sample_batches = sample_batches,
            sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
            runtime_attr_rdtest=runtime_attr_rdtest
        }
    
    call igv_individual.IGV_all_samples as igv_plots {
        input:
            prefix = prefix,
            varfile = varfile,
            Fasta = Fasta,
            Fasta_dict = Fasta_dict,
            Fasta_idx = Fasta_idx,
            samples = sample_list,
            crams = cram_list,
            crams_idx = crai_list,
            sv_base_mini_docker = sv_base_mini_docker,
            igv_docker = igv_docker,
            runtime_attr_override=runtime_attr_override
        }

    call concatinate_plots{
        input:
            rd_plots = RdTest.Plots,
            igv_plots = igv_plots.tar_gz_pe,
            prefix = prefix,
            varfile = varfile,
            pedfile = pedfile,
            igv_docker = igv_docker,
            runtime_attr_concatinate = runtime_attr_concatinate
    }

    output{
        File concatinated_plots = concatinate_plots.plots
    }
}

task concatinate_plots{
    input{
        File rd_plots
        File igv_plots
        String prefix
        File varfile
        File pedfile
        String igv_docker
        RuntimeAttr? runtime_attr_concatinate
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_concatinate, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: igv_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    command <<<
        set -eu -o pipefail    

        tar -zxf ~{rd_plots}
        tar -zxf ~{igv_plots}
        mkdir ~{prefix}_igv_rdtest_plots
        echo 'test'
        python3 /src/MakeRDtest.py \
            ~{varfile} \
            ~{pedfile} \
            ~{prefix} \
            10000000 \
            ~{prefix}_igv_plots \
            ~{prefix}_rd_plots \
            ~{prefix}_igv_rdtest_plots
        tar -czf ~{prefix}_igv_rdtest_plots.tar.gz ~{prefix}_igv_rdtest_plots
    >>>

    output{
        File plots = "~{prefix}_igv_rdtest_plots.tar.gz"
    }


}

