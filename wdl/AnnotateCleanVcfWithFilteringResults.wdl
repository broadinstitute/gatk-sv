version 1.0

#script to annotate cleanvcf with the results from downstream filterings

import "Structs.wdl"

workflow IntegrateFilteringResultsToCleanVcf {
    input{
        File clean_vcf
        File clean_vcf_idx
        File minGQ_vcf
        File minGQ_vcf_idx
        File outlier_removal_vcf
        File outlier_removal_vcf_idx
        File batch_effect_vcf
        File batch_effect_vcf_idx
        File final_recali_vcf
        File final_recali_vcf_idx
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_anno_FinalRecali
        RuntimeAttr? runtime_attr_anno_BatchEffect
        RuntimeAttr? runtime_attr_anno_OutlierFilter
        RuntimeAttr? runtime_attr_anno_minGQ
        RuntimeAttr? runtime_attr_ConcatVcfs
        }

    call AnnotateFinalRecali{
        input:
            vcf = batch_effect_vcf,
            filtered_vcf = final_recali_vcf,
            sv_pipeline_docker =  sv_pipeline_docker,
            runtime_attr_override = runtime_attr_anno_FinalRecali
    }

    call AnnotateBatchEffect{
        input:
            vcf = outlier_removal_vcf,
            filtered_vcf = AnnotateFinalRecali.output_vcf,
            sv_pipeline_docker =  sv_pipeline_docker,
            runtime_attr_override = runtime_attr_anno_BatchEffect
    }

    call AnnotateOutlierFilter{
        input:
            vcf = minGQ_vcf,
            filtered_vcf = AnnotateBatchEffect.output_vcf,
            sv_pipeline_docker =  sv_pipeline_docker,
            runtime_attr_override = runtime_attr_anno_OutlierFilter
    }

    call AnnotateMinGQ{
        input:
            vcf = clean_vcf,
            filtered_vcf = AnnotateOutlierFilter.output_vcf,
            prefix = prefix,
            sv_pipeline_docker =  sv_pipeline_docker,    
            runtime_attr_override = runtime_attr_anno_minGQ
    }

    output{
        File annotated_vcf = AnnotateMinGQ.output_vcf
        File annotated_vcf_idx = AnnotateMinGQ.output_vcf_idx
    }
}

task AnnotateFinalRecali{
    input{
        File vcf
        File filtered_vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euo pipefail

        zcat ~{filtered_vcf} |fgrep -v "#" |cut -f3,7 > SVID_filter

        python <<CODE
        import os

        fin=open("SVID_filter")
        svid_hash = {}
        for line in fin:
            pin=line.strip().split()
            svid_hash[pin[0]]=pin[1]
        fin.close()

        fin=os.popen(r'''zcat %s'''%("~{filtered_vcf}"))
        header = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': 
                header.append(pin)
            else:
                break
        fin.close()

        header.append(['##FILTER=<ID=FAIL_FINAL_RECALIBRATION,Description="SV failed final cleanup and quality recalibration">'])


        fo=open("annotated.vcf",'w')
        for i in header:
            print(' '.join(i), file=fo)

        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': continue
            elif pin[0][0]=='#': print('\t'.join(pin), file=fo)
            else:
                if pin[2] in svid_hash.keys():
                    pin[6]=svid_hash[pin[2]]
                else:
                    pin[6]='FAIL_FINAL_RECALIBRATION'
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip annotated.vcf
        tabix annotated.vcf.gz
    >>>

    output{
        File output_vcf = "annotated.vcf.gz"
        File output_vcf_idx = "annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task AnnotateBatchEffect{
    input{
        File vcf
        File filtered_vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euo pipefail

        zcat ~{filtered_vcf} | grep -v "#" | cut -f3,7 > SVID_filter

        python <<CODE
        import os

        fin=open("SVID_filter")
        svid_hash = {}
        for line in fin:
            pin=line.strip().split()
            svid_hash[pin[0]]=pin[1]
        fin.close()

        fin=os.popen(r'''zcat %s'''%("~{filtered_vcf}"))
        header = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': 
                header.append(pin)
            else:
                break
        fin.close()

        header.append(['##FILTER=<ID=FAIL_BATCH_EFFECT,Description="SV failed batch effect check">'])

        fo=open("annotated.vcf",'w')
        for i in header:
            print(' '.join(i), file=fo)

        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': continue
            elif pin[0][0]=='#': print('\t'.join(pin), file=fo)
            else:
                if pin[2] in svid_hash.keys():
                    pin[6]=svid_hash[pin[2]]
                else:
                    pin[6]='FAIL_BATCH_EFFECT'
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip annotated.vcf
        tabix annotated.vcf.gz
    >>>

    output{
        File output_vcf = "annotated.vcf.gz"
        File output_vcf_idx = "annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task AnnotateOutlierFilter{
    input{
        File vcf
        File filtered_vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euo pipefail

        zcat ~{filtered_vcf} | grep -v "#" | cut -f3,7 > SVID_filter

        python <<CODE
        import os

        fin=open("SVID_filter")
        svid_hash = {}
        for line in fin:
            pin=line.strip().split()
            svid_hash[pin[0]]=pin[1]
        fin.close()

        fin=os.popen(r'''zcat %s'''%("~{filtered_vcf}"))
        header = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': 
                header.append(pin)
            else:
                break
        fin.close()

        header.append(['##FILTER=<ID=FAIL_OUTLIER_REMOVAL,Description="SV failed outlier removal">'])

        fo=open("annotated.vcf",'w')
        for i in header:
            print(' '.join(i), file=fo)

        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': continue
            elif pin[0][0]=='#': print('\t'.join(pin), file=fo)
            else:
                if pin[2] in svid_hash.keys():
                    pin[6]=svid_hash[pin[2]]
                else:
                    pin[6]='FAIL_OUTLIER_REMOVAL'
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip annotated.vcf
        tabix annotated.vcf.gz
    >>>

    output{
        File output_vcf = "annotated.vcf.gz"
        File output_vcf_idx = "annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task AnnotateMinGQ{
    input{
        File vcf
        File filtered_vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euo pipefail

        zcat ~{filtered_vcf} | grep -v "#" | cut -f3,7 > SVID_filter

        python <<CODE
        import os

        fin=open("SVID_filter")
        svid_hash = {}
        for line in fin:
            pin=line.strip().split()
            svid_hash[pin[0]]=pin[1]
        fin.close()

        fin=os.popen(r'''zcat %s'''%("~{filtered_vcf}"))
        header = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': 
                header.append(pin)
            else: 
                break
        fin.close()

        header.append(['##FILTER=<ID=FAIL_minGQ,Description="SV failed minGQ">'])

        fo=open("annotated.vcf",'w')
        for i in header:
            print(' '.join(i), file=fo)

        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': continue
            elif pin[0][0]=='#': print('\t'.join(pin), file=fo)
            else:
                if pin[2] in svid_hash.keys():
                    pin[6]=svid_hash[pin[2]]
                else:
                    pin[6]='FAIL_minGQ'
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip annotated.vcf
        mv annotated.vcf.gz ~{prefix}.filter_annotated.vcf.gz
        tabix ~{prefix}.filter_annotated.vcf.gz
    >>>

    output{
        File output_vcf = "~{prefix}.filter_annotated.vcf.gz"
        File output_vcf_idx = "~{prefix}.filter_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task AnnotateWithFilterResults{
    input{
        File vcf
        File filtered_vcf
        String prefix
        String new_header
        String new_filter
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euo pipefail

        zcat ~{filtered_vcf} | grep -v "#" | cut -f3,7 > SVID_filter

        python <<CODE
        import os

        fin=open("SVID_filter")
        svid_hash = {}
        for line in fin:
            pin=line.strip().split()
            svid_hash[pin[0]]=pin[1]
        fin.close()

        fin=os.popen(r'''zcat %s'''%("~{filtered_vcf}"))
        header = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': 
                header.append(pin)
        fin.close()

        header.append('~{new_header}')

        fo=open("annotated.vcf",'w')
        for i in header:
            print(' '.join(i), file=fo)

        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##': continue
            elif pin[0][0]=='#': print('\t'.join(pin), file=fo)
            else:
                if pin[2] in svid_hash.keys():
                    pin[6]=svid_hash[pin[2]]
                else:
                    pin[6]='~{new_filter}'
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip annotated.vcf
        tabix annotated.vcf.gz
    >>>

    output{
        File output_vcf = "~{prefix}.filter_annotated.vcf.gz"
        File output_vcf_idx = "~{prefix}.filter_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
