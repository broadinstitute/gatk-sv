version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks


workflow MendelianInheritanceStats {

    input {
        Array[File] vcfs
        File ped
        File inheritance_table   # fa mo pb category
        String output_prefix 
        String sv_pipeline_base_docker

        File anno_script_bash
        File anno_script_helper_R
        File repeat_mask
        File simple_repeats 
        File segmental_duplicates

        RuntimeAttr? runtime_attr_identify_families
        RuntimeAttr? runtime_attr_normalize_vcf
        RuntimeAttr? runtime_attr_extract_fam_non_ref
        RuntimeAttr? runtime_attr_classify_variants
        RuntimeAttr? runtime_attr_median_counts
        RuntimeAttr? runtime_attr_merge_counts
    }

    scatter (vcf in vcfs) {

        call IdentifyFamilies {
            input:
                ped = ped,
                vcf = vcf,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_identify_families
        }

        call NormalizeVCF {
            input:
                vcf = vcf,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_normalize_vcf
        }

        call ExtractFamilyNonRef {
            input:
                vcf = NormalizeVCF.out_vcf,
                families = IdentifyFamilies.family_table,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_extract_fam_non_ref
        }


       call LongReadGenotypeTasks.AnnotateGenomicContext{
          input:
            variant_sites = ExtractFamilyNonRef.family_vcf,
            anno_script_bash = anno_script_bash,
            anno_script_Rscript = anno_script_helper_R,
            repeat_mask = repeat_mask,
            simple_repeats = simple_repeats,
            segmental_duplicates = segmental_duplicates,
            docker_image = sv_pipeline_base_docker
        }

 
        call ClassifyVariants {
            input:
                vcf = ExtractFamilyNonRef.family_vcf,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_classify_variants
        }

        call MendelianCounts {
            input:
                vcf = ClassifyVariants.classified_vcf,
                inheritance_table = inheritance_table,
                families = IdentifyFamilies.family_table,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_median_counts
        }
    }

    call MergeCounts {
        input:
            count_tables = MendelianCounts.count_table,
            prefix = output_prefix,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_merge_counts
    }

    output {
        File final_counts = MergeCounts.counts
        File final_proportions = MergeCounts.proportions
    }
}

task IdentifyFamilies {

    input {
        File ped
        File vcf
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        bcftools query -l ~{vcf} > samples.txt

        python <<EOF
        
        import pandas as pd

        ped = pd.read_csv("~{ped}", sep="\t", header=None)
        ped.columns=["fam","id","father","mother","sex","pheno"]

        samples=set(open("samples.txt").read().split())

        out=open("families.txt","w")

        for fam,df in ped.groupby("fam"):

            for _,row in df.iterrows():

                if row.father!="0" and row.mother!="0":

                    if row.id in samples and row.father in samples and row.mother in samples:
                        out.write(f"{fam}\t{row.father}\t{row.mother}\t{row.id}\n")

        out.close()
        EOF
    >>>

    output {
        File family_table = "families.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task NormalizeVCF {

    input {
        File vcf
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        bcftools norm -m -any ~{vcf} -Oz -o norm.vcf.gz
        bcftools index norm.vcf.gz
    >>>

    output {
        File out_vcf = "norm.vcf.gz"
        File out_index = "norm.vcf.gz.csi"
    }


    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractFamilyNonRef {

    input {
        File vcf
        File families
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f2- ~{families} | sed -e 's/\t/\n/g' > sample_list

        bcftools view -S sample_list -c 1 ~{vcf} -Oz -o family_nonref.vcf.gz

        bcftools index family_nonref.vcf.gz
    >>>

    output {
        File family_vcf = "family_nonref.vcf.gz"
        File family_vcf_index = "family_nonref.vcf.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ClassifyVariants {

    input {
        File vcf
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python <<EOF
        import pysam

        vcf=pysam.VariantFile("~{vcf}")

        header=vcf.header
        header.add_line('##INFO=<ID=VARCLASS,Number=1,Type=String,Description="variant class">')

        out=pysam.VariantFile("classified.vcf","w",header=header)

        for rec in vcf:

            ref=len(rec.ref)
            alt=len(rec.alts[0])

            if ref==1 and alt==1:
                cls="SNV"

            elif alt>ref:
                l=alt-ref
                if l<=29: cls="INS_1_29"
                elif l<=49: cls="INS_30_49"
                else: cls="INS_50P"

            else:
                l=ref-alt
                if l<=29: cls="DEL_1_29"
                elif l<=49: cls="DEL_30_49"
                else: cls="DEL_50P"

            rec.info["VARCLASS"]=cls
            out.write(rec)

        out.close()
        EOF

        bgzip classified.vcf
        tabix -p vcf classified.vcf.gz
    >>>

    output {
        File classified_vcf = "classified.vcf.gz"
        File classified_index = "classified.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MendelianCounts {

    input {
        File vcf
        File inheritance_table
        File families
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python <<EOF
        import pysam
        import pandas as pd
        from collections import defaultdict

        def norm(gt):
            if gt is None:
                return "./."
            gt = str(gt).split(":")[0]   # remove FORMAT fields
            gt = gt.replace("|", "/")    # remove phasing
            alleles = gt.split("/")
            if len(alleles) != 2:
                return None
            a, b = sorted(alleles)
            return f"{a}/{b}"

        def norm_individual_gt(gt):
            if gt is None:
                return "./."
            return "/".join(map(str,sorted(gt)))



        families=[l.strip().split() for l in open("~{families}")]

        rules = pd.read_csv("~{inheritance_table}", sep="\t")

        rules_dict = {
            (norm(r.fa), norm(r.mo), norm(r.pb)): r.category
            for _, r in rules.iterrows()
        }


        counts=defaultdict(int)

        vcf=pysam.VariantFile("~{vcf}")
        for rec in vcf:
            vclass=rec.info["VARCLASS"]
            for fam,fa,mo,pb in families:
                f = norm_individual_gt(rec.samples[fa]["GT"])
                m = norm_individual_gt(rec.samples[mo]["GT"])
                p = norm_individual_gt(rec.samples[pb]["GT"])
                key=(f,m,p)
                if key in rules_dict:
                    cat=rules_dict[key]
                    counts[(vclass,cat,fam)]+=1

        with open("counts.tsv","w") as out:
            out.write("variant_class\tcategory\tfamily\tcount\n")
            for (v,c,f),n in counts.items():
                out.write(f"{v}\t{c}\t{f}\t{n}\n")


        EOF
    >>>

    output {
        File count_table = "counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeCounts {

    input {
        Array[File] count_tables
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        cat ~{sep=" " count_tables} > all_counts.tmp

        python <<EOF

        import pandas as pd
        import sys

        # read table
        df = pd.read_csv("all_counts.tmp", sep="\t")

        # convert count to numeric; non-numeric values become NaN
        df["count"] = pd.to_numeric(df["count"], errors="coerce")

        # keep only rows where count is a valid number
        df = df[df["count"].notna()]


        for fam, fam_df in df.groupby("family"):

            # sum counts per (variant_class, category)
            summary = (
                fam_df.groupby(["variant_class", "category"], as_index=False)["count"]
                .sum()
            )

            # calculate proportion within variant_class
            summary["total"] = summary.groupby("variant_class")["count"].transform("sum")
            summary["prop"] = (summary["count"] / summary["total"]).round(3)

            # pivot tables
            count_table = summary.pivot(
                index="variant_class",
                columns="category",
                values="count"
            ).fillna(0).astype(int)

            prop_table = summary.pivot(
                index="variant_class",
                columns="category",
                values="prop"
            ).fillna(0)

            # write outputs
            count_table.to_csv(f"~{prefix}.{fam}.count.tsv", sep="\t")
            prop_table.to_csv(f"~{prefix}.{fam}.prop.tsv", sep="\t")

        EOF
    >>>

    output {
        File counts = "~{prefix}.count.tsv"
        File proportions = "~{prefix}.prop.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



