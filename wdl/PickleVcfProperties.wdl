version 1.0

import "Utils.wdl" as Utils

workflow PickleVcfProperties {
    input {
        File vcf
        File vcf_index
        Array[String] wanted_properties
        String samtools_cloud_docker
        String sv_utils_docker
    }

    call Utils.GetVcfSize {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call GetNeededMemGB {
        input:
            wanted_properties=wanted_properties,
            num_entries=GetVcfSize.num_entries
    }

    call ReadAndPickleProperties {
        input:
          vcf=vcf,
          wanted_properties=wanted_properties,
          mem_gb=GetNeededMemGB.pickle_task_mem_gb,
          sv_utils_docker=sv_utils_docker
    }

    output {
        Int num_records = GetVcfSize.num_records
        Int num_samples = GetVcfSize.num_samples
        Float num_entries = GetVcfSize.num_entries
        File pickled_properties = ReadAndPickleProperties.pickled_properties
    }
}

task GetNeededMemGB {
    input {
        Array[String] wanted_properties
        Float num_entries
    }

    Float mem_gb = "1.5"
    Int disk_gb = 10

    runtime {
        docker: "python:slim"
        cpu: 1
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
    python - <<'CODE'
properties=["~{sep='", "' wanted_properties}"]
entries_cost_map={"svtype": 0.0, "svlen": 0.0, "gt": 1.4e-8, "gq": 5e-8, "sl": 5e-8}
mem_gb = ~{num_entries} * sum(entries_cost_map.get(property, 5e-8) for property in properties)
print(mem_gb)
CODE
    >>>

    output {
        Float pickle_task_mem_gb = read_float(stdout())
    }
}

task ReadAndPickleProperties {
    input {
        File vcf
        Array[String] wanted_properties
        String sv_utils_docker
        Float mem_gb
    }

    # High disk size for large throughput. A large proportion of run time is loading data from huge VCFs. Disk is cheap.
    Int disk_gb = round(
        1000 + size([vcf], "GiB")
    )

    # create output filename:
    #   strip .vcf.gz from end (if present)
    #   add wanted properties (separated by underscores)
    #   and add ".pickle.bz2"
    #   (sep doesn't work the way it's supposed to outside of command block, so move part of this logic to command)
    String vcf_basename = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "")

    runtime {
        docker: sv_utils_docker
        cpu: 1
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<

    python - <<'CODE'
from sv_utils import benchmark_variant_filter

benchmark_variant_filter.load_benchmark_properties_from_vcf_properties(
  "~{vcf}", wanted_properties=["~{sep='", "' wanted_properties}"]
).to_pickle("~{vcf_basename}_~{sep='_' wanted_properties}.pickle")
CODE
    >>>

    output {
        File pickled_properties = select_first(glob("*.pickle"))
    }
}