version 1.0

import "Structs.wdl"

workflow RdTestVisualization{
    input{
        String prefix
        Array[File] medianfile
        File pedfile
        File sample_batches
        File batch_bincov
        File bed
        String flags
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_rdtest
    }
        call rdtest{
            input:
                bed=bed,
                medianfile=medianfile,
                pedfile=pedfile,
                sample_batches=sample_batches,
                batch_bincov=batch_bincov,
                prefix=prefix,
                flags=flags,
                sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
        output{
            File Plots = rdtest.plots
        }
}


# Run rdtest
task rdtest {
    input{
        File bed
        File sample_batches # samples, batches
        File batch_bincov # batch, bincov
        Array[File] medianfile
        File pedfile
        String prefix
        String sv_pipeline_rdtest_docker
        String flags
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    command <<<
        set -ex
        awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$5}' ~{bed} |egrep "DEL|DUP" | sort -k1,1 -k2,2n> test.bed
        cut -f5 test.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        fgrep -wf samples.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
        fgrep -f existing_batches.txt ~{batch_bincov} > bincovlist.txt
        paste ~{sep=" " medianfile} > medianfile.txt

        i=0
        bedtools merge -i test.bed > test.merged.bed
        while read batch bincov
        do
            let "i=$i+1"
            if [ $i -gt 1 ]
            then
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R test.merged.bed|cut -f4->covfile.$i.bed 
            else
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R test.merged.bed>covfile.$i.bed 

        fi
        done<bincovlist.txt

        paste covfile.*.bed |tr ' ' '\t' |bgzip >allcovfile.bed.gz 
        tabix allcovfile.bed.gz
        rm covfile.*.bed
        zcat allcovfile.bed.gz |head -n 1|cut -f 4-|tr '\t' '\n'>samples.txt
        Rscript /opt/RdTest/RdTest.R \
            -b test.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -f ~{pedfile} \
            -p TRUE \
            -w samples.txt \
            -s 10000000 \
            ~{flags}
        mkdir ~{prefix}_rd_plots
        mv *jpg ~{prefix}_rd_plots
        tar -czvf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
    >>>
    
    output {
        File plots = "~{prefix}_rd_plots.tar.gz"
    }
    
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_rdtest_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
