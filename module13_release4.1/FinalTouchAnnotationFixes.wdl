#this script is developed to fix false postiive inversions in gnomad v3 vcf that are >1Mb in size but have PE_GT=0
#INVs with no passing samples have filter column revised from "PASS" to "UNRESOLVED"

#developing workdir on erisone: /data/talkowski/xuefang/data/gnomad_V3/module08/step9_sm_depyh_only_dup_fix

version 1.0

import "Structs.wdl"

workflow FinalTouchAnnotationFixes {

    input {
        Array[File] vcfs
        Array[File] vcf_idxes
        File SVID_reference_artifacts #a list of SVID that are likely reference artifacts. over 99% of samples carrying these sites are hom-alts

        File? subsets_af_header_lines
        Array[File?] subsets_af_info_list
        Boolean fix_mCNV_formats = false

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_recognize_ref_arti
        RuntimeAttr? runtime_attr_revise_cpx_func_anno
    }

    scatter (i in range(length(vcfs))){

        call RecognizeReferenceArtifacts{
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxes[i],
                SVID_reference_artifacts = SVID_reference_artifacts,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_recognize_ref_arti
        }

        call ReviseCpxFunctionalAnnotations{
            input:
                vcf = RecognizeReferenceArtifacts.reference_artifact_annotated_vcf,
                vcf_idx = RecognizeReferenceArtifacts.reference_artifact_annotated_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_revise_cpx_func_anno
        }

        call ReviseMcnvInfo{
            input:
                vcf = ReviseCpxFunctionalAnnotations.cpx_anno_revised_vcf,
                vcf_idx = ReviseCpxFunctionalAnnotations.cpx_anno_revised_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call ReviseCpxIntervals{
            input:
                vcf = ReviseMcnvInfo.mCNV_filter_cleaned_vcf,
                vcf_idx = ReviseMcnvInfo.mCNV_filter_cleaned_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        if(fix_mCNV_formats){
            call ReviseMcnvFormat{
                input:
                    vcf = ReviseCpxIntervals.cpx_interval_revised_vcf,
                    vcf_idx = ReviseCpxIntervals.cpx_interval_revised_vcf_idx,
                    sv_pipeline_base_docker = sv_pipeline_base_docker
            }
        }

        File fixed_vcf = select_first([ReviseMcnvFormat.mCNV_format_cleaned_vcf, ReviseCpxIntervals.cpx_interval_revised_vcf])
        File fixed_vcf_idx = select_first([ReviseMcnvFormat.mCNV_format_cleaned_vcf_idx, ReviseCpxIntervals.cpx_interval_revised_vcf_idx])

        call ExtractVcfSites as extract_vcf_sites_1{
            input:
                vcf = fixed_vcf,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call Vcf2Bed as vcf_to_bed{
            input:
                vcf = fixed_vcf,
                vcf_idx = fixed_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        if (defined(subsets_af_header_lines)){
            call AddAfOfSubsets{
                input:
                    vcf = fixed_vcf,
                    vcf_idx = fixed_vcf_idx,
                    subsets_af_header_lines = subsets_af_header_lines,
                    subsets_af_info = subsets_af_info_list[i],
                    sv_pipeline_base_docker = sv_pipeline_base_docker
            }

            call ExtractVcfSites as extract_vcf_sites_2{
                input:
                    vcf = AddAfOfSubsets.subset_af_annotated_vcf,
                    sv_pipeline_base_docker = sv_pipeline_base_docker
            }
        }
    }

    output{
        Array[File] revised_vcf          = fixed_vcf
        Array[File] revised_vcf_idx      = fixed_vcf_idx
        Array[File] revised_sites        = extract_vcf_sites_1.vcf_sites
        Array[File] revised_sites_idx    = extract_vcf_sites_1.vcf_sites_idx
        Array[File] revised_bed          = vcf_to_bed.zipped_bed
        Array[File?] subset_af_vcf       = AddAfOfSubsets.subset_af_annotated_vcf
        Array[File?] subset_af_vcf_idx   = AddAfOfSubsets.subset_af_annotated_vcf_idx
        Array[File?] subset_af_sites     = extract_vcf_sites_2.vcf_sites
        Array[File?] subset_af_sites_idx = extract_vcf_sites_2.vcf_sites_idx
    }
}

task RecognizeReferenceArtifacts{
    input{
        File vcf
        File vcf_idx
        File SVID_reference_artifacts

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        python <<CODE

        import os
        import pysam

        #readin the SVID of reference artifacts SVs
        fin = open("~{SVID_reference_artifacts}")
        SVID_list = []
        for line in fin:
            pin=line.strip().split()
            SVID_list += pin
        fin.close()

        fin = pysam.VariantFile("~{vcf}")
        header = fin.header
        header.add_line('##FILTER=<ID=REFERENCE_ARTIFACT,Description="Likely reference artifact sites that are homozygous alternative in over 99% of the samples">')
        header.info.remove_header('LIKELY_REFERENCE_ARTIFACT')

        fo = pysam.VariantFile("~{prefix}.ref_arti_annotated.vcf.gz", 'w', header = header)
        for record in fin:
            if record.id in SVID_list:
                if 'PASS' in record.filter.keys():
                    record.filter.add('REFERENCE_ARTIFACT')
            if 'LIKELY_REFERENCE_ARTIFACT' in record.info.keys():
                del record.info['LIKELY_REFERENCE_ARTIFACT']
            fo.write(record)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.ref_arti_annotated.vcf.gz"

    >>>

    output{
        File reference_artifact_annotated_vcf = "~{prefix}.ref_arti_annotated.vcf.gz"
        File reference_artifact_annotated_vcf_idx = "~{prefix}.ref_arti_annotated.vcf.gz.tbi"
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

task ReviseCpxFunctionalAnnotations{
    input{
        File vcf
        File vcf_idx
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        python <<CODE

        import os
        import pysam


        fin = pysam.VariantFile("~{vcf}")
        header = fin.header
        header.add_line('##INFO=<ID=PREDICTED_COMPLEX_CODING_DISRUPTION,Number=.,Type=String,Description="Complex SVs that interrupt genes">')
        fo = pysam.VariantFile("~{prefix}.cpx_anno_revised.vcf.gz", 'w', header = header)

        for record in fin:
            if record.info['SVTYPE'] == "CPX":
                interrupted_genes = ()
                for anno in record.info.keys():
                    if "PREDICTED_" in anno and not anno=='PREDICTED_INTERGENIC':
                        interrupted_genes+=record.info[anno]
                        del record.info[anno]
                if len(interrupted_genes)>0:
                    record.info['PREDICTED_COMPLEX_CODING_DISRUPTION'] = interrupted_genes
            fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.cpx_anno_revised.vcf.gz"

    >>>

    output{
        File cpx_anno_revised_vcf = "~{prefix}.cpx_anno_revised.vcf.gz"
        File cpx_anno_revised_vcf_idx = "~{prefix}.cpx_anno_revised.vcf.gz.tbi"
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

task ReviseMcnvInfo{
    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        python <<CODE

        import os
        import pysam

        #readin the SVID of reference artifacts SVs

        fin = pysam.VariantFile("~{vcf}")
        header = fin.header
        header.filters.remove_header('MULTIALLELIC')

        fo = pysam.VariantFile("~{prefix}.mCNV_filter_cleaned.vcf.gz", 'w', header = header)
        for record in fin:
            if 'MULTIALLELIC' in record.filter.keys():
                del(record.filter['MULTIALLELIC'])
            fo.write(record)
        
        fin.close()
        fo.close()


        CODE

        tabix -p vcf "~{prefix}.mCNV_filter_cleaned.vcf.gz"

    >>>

    output{
        File mCNV_filter_cleaned_vcf = "~{prefix}.mCNV_filter_cleaned.vcf.gz"
        File mCNV_filter_cleaned_vcf_idx = "~{prefix}.mCNV_filter_cleaned.vcf.gz.tbi"
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

task ReviseMcnvFormat{
    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        

        bcftools view -h ~{vcf} | bgzip > ~{prefix}.non_mCNV.vcf.gz
        bcftools view -H ~{vcf} | awk '{if ($5!="<CNV>") print}' | bgzip >> ~{prefix}.non_mCNV.vcf.gz
        bcftools view -h ~{vcf} | bgzip > ~{prefix}.mCNV.vcf.gz
        bcftools view -H ~{vcf} | awk '{if ($5=="<CNV>") print}' | bgzip >> ~{prefix}.mCNV.vcf.gz


        python <<CODE

        import os
        import gzip

        fin = os.popen(r'''zcat %s'''%("~{prefix}.mCNV.vcf.gz"))
        fo = open("~{prefix}.mCNV_format_cleaned.vcf", 'w')
        for line in fin:
            pin = line.strip().split()
            if pin[0][:2] == "##":
                print(' '.join(pin), file = fo)
            elif pin[0][0] == "#":
                print('\t'.join(pin), file = fo)
            else:
                if pin[4] == "<CNV>" and not "CNQ" in pin[8]:
                    format_new = pin[:8]
                    format_key = [pin[8].split(':')[0]] + ['CN','CNQ'] + pin[8].split(':')[1:]
                    gq_pos = format_key.index('GQ')
                    del format_key[gq_pos]
                    format_new.append(':'.join(format_key))
                    for i in pin[9:]:
                        sample_new = [i.split(':')[0], i.split(':')[pin[8].split(':').index('RD_CN')],i.split(':')[pin[8].split(':').index('RD_GQ')]] + i.split(':')[1:]
                        del sample_new[gq_pos]
                        format_new.append(':'.join(sample_new))
                    print('\t'.join(format_new), file = fo)
                else:
                    print('\t'.join(pin), file = fo)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.non_mCNV.vcf.gz" 
        bgzip "~{prefix}.mCNV_format_cleaned.vcf"
        tabix -p vcf "~{prefix}.mCNV_format_cleaned.vcf.gz"
        echo "~{prefix}.non_mCNV.vcf.gz" > vcf_list.tsv
        echo "~{prefix}.mCNV_format_cleaned.vcf.gz" >> vcf_list.tsv
        bcftools concat --allow-overlaps --file-list vcf_list.tsv -O z -o ~{prefix}.mCNV_format_organized.vcf.gz
        tabix -p vcf "~{prefix}.mCNV_format_organized.vcf.gz"

    >>>

    output{
        File mCNV_format_cleaned_vcf = "~{prefix}.mCNV_format_organized.vcf.gz"
        File mCNV_format_cleaned_vcf_idx = "~{prefix}.mCNV_format_organized.vcf.gz.tbi"
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


task Vcf2Bed{
    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 2),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        svtk vcf2bed -i ALL --include-filters ~{vcf} ~{prefix}.bed
        bgzip "~{prefix}.bed"

    >>>

    output{
        File zipped_bed = "~{prefix}.bed.gz"
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


#this function is built to re-format the INS-iDEL CPX events
#in gnomadv3, the INS_iDEL has CPX_interal=DEL_chrN:pos1-end1,INS_chrN:pos2-end2
#with this fix, the END2 and POS2 will be removed from CPX events;   INS_iDEL will have SOURCE=INS_chrN:pos2-end2, and CPX_interal=DEL_chrN:pos1-end1 
task ReviseCpxIntervals{

    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        python <<CODE

        import os
        import pysam

        #readin the SVID of reference artifacts SVs

        fin = pysam.VariantFile("~{vcf}")
        header = fin.header

        fo = pysam.VariantFile("~{prefix}.cpx_interval_revised.vcf.gz", 'w', header = header)
        for record in fin:
            if record.info['SVTYPE']=="CPX":
                if "POS2" in record.info.keys():
                    del record.info['POS2']
                if "END2" in record.info.keys():
                    del record.info['END2']
                if record.info['CPX_TYPE'] == "INS_iDEL":
                    DEL_section = [i for i in record.info['CPX_INTERVALS'] if i.split('_')[0]=='DEL'][0]
                    INS_section = [i for i in record.info['CPX_INTERVALS'] if i.split('_')[0]=='INS'][0]
                    record.info['SOURCE'] = INS_section
                    record.info['CPX_INTERVALS'] = (DEL_section)
            fo.write(record)
        
        fin.close()
        fo.close()


        CODE

        tabix -p vcf "~{prefix}.cpx_interval_revised.vcf.gz"

    >>>

    output{
        File cpx_interval_revised_vcf = "~{prefix}.cpx_interval_revised.vcf.gz"
        File cpx_interval_revised_vcf_idx = "~{prefix}.cpx_interval_revised.vcf.gz.tbi"
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

#add the af info of subsets to the vcf
task AddAfOfSubsets{

    input{
        File vcf
        File vcf_idx
        File? subsets_af_header_lines
        File? subsets_af_info

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: ceil(20.0 + size(vcf, "GiB") * 2), 
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        python <<CODE

        import os
        import pysam

        def readin_af_info_of_subsets(af_subsets):
            af_subsets_info = {}
            fin = os.popen(r'''zcat %s'''%(af_subsets))
            for line in fin:
                pin=line.strip().split()
                af_subsets_info[pin[0]] = pin[1:]
            fin.close()
            return af_subsets_info

        def readin_af_header(header_lines):
            AF_subsets_INFO = []
            fin = open(header_lines)
            for line in fin:
                pin=line.strip().split('\t')
                AF_subsets_INFO+=(pin)

            fin.close()
            return AF_subsets_INFO

        af_subsets_info = readin_af_info_of_subsets("~{subsets_af_info}")
        AF_subsets_INFO = readin_af_header("~{subsets_af_header_lines}")

        fin=pysam.VariantFile("~{vcf}")
        for line in AF_subsets_INFO:
            fin.header.add_line(line)

        fo = pysam.VariantFile("~{prefix}.with_subset_af.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if not record.id in af_subsets_info.keys():
                print(record.id)
                break
            for i in range(len(af_subsets_info['name'])):
                if af_subsets_info[record.id][i]=="NA": continue
                if "," in af_subsets_info[record.id][i]:
                    if "_FREQ" in af_subsets_info['name'][i]:
                        record.info[af_subsets_info['name'][i]] = tuple([float(j) for j in af_subsets_info[record.id][i].split(',')])
                    else:
                        record.info[af_subsets_info['name'][i]] = tuple([int(j) for j in af_subsets_info[record.id][i].split(',')])
                elif af_subsets_info['name'][i].split('_')[-1] == 'AC':
                    record.info[af_subsets_info['name'][i]] = (int(af_subsets_info[record.id][i]),)
                elif af_subsets_info['name'][i].split('_')[-1] == 'AF':
                    record.info[af_subsets_info['name'][i]] = (float(af_subsets_info[record.id][i]),)
                elif "_FREQ" in af_subsets_info['name'][i]:
                    record.info[af_subsets_info['name'][i]] = float(af_subsets_info[record.id][i])
                else: 
                    record.info[af_subsets_info['name'][i]] = int(af_subsets_info[record.id][i])
            fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.with_subset_af.vcf.gz"
    >>>

    output{
        File subset_af_annotated_vcf = "~{prefix}.with_subset_af.vcf.gz"
        File subset_af_annotated_vcf_idx = "~{prefix}.with_subset_af.vcf.gz.tbi"
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

task ExtractVcfSites{
    input{
        File vcf
        File? vcf_idx
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(20.0 + size(vcf, "GiB") * 1.5),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        
        zcat ~{vcf} | cut -f1-8 | bgzip > "~{prefix}.sites.gz"
        tabix -p vcf "~{prefix}.sites.gz"
    >>>

    output{
        File vcf_sites = "~{prefix}.sites.gz"
        File vcf_sites_idx = "~{prefix}.sites.gz.tbi"
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



