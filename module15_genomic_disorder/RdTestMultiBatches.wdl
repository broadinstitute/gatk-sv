#script to generate rdplots across multiple batches
version 1.0

import "Structs.wdl"

workflow RdTestMultiBatches {

	input{
		File rdtest_bed
    File? whitelist
		String prefix
		File sample_batch_list
		File Rdtest_V2_script
		String sv_pipeline_base_docker
	}


	call RemoteTabixRdMetrics{
		input:
			rdtest_bed = rdtest_bed,
			sample_batch_list = sample_batch_list,
			sv_pipeline_base_docker = sv_pipeline_base_docker
	}


	call CollectRdMedian{
		input:
			rdtest_bed = rdtest_bed,
			sample_batch_list = sample_batch_list,
			sv_pipeline_base_docker = sv_pipeline_base_docker
	}

	call RdTest{
		input:
			rdtest_bed = rdtest_bed,
 			prefix = prefix,
      whitelist = whitelist,
			rd_metrics_folder = RemoteTabixRdMetrics.rd_metrics_folder,
			rd_median_file = CollectRdMedian.rd_median_file,
			Rdtest_V2_script  = Rdtest_V2_script,
			sv_pipeline_base_docker = sv_pipeline_base_docker

	}

  output{
    File rd_plot_tarball = RdTest.rd_plots
  }
}

task RdTest {
  input {
    File rdtest_bed
    File? whitelist
    String prefix
    File rd_metrics_folder
    File rd_median_file
    File Rdtest_V2_script
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: ceil(3 + size(rd_metrics_folder, "GB") * 5),
    disk_gb: ceil(10 + size(rd_metrics_folder, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  String whitelist_arg = if defined(whitelist) then "-w whitelist.local" else ""
  command <<<

    set -euo pipefail

    if [[ "~{defined(whitelist)}" == "true" ]]; then
      gsutil cp "~{whitelist}" whitelist.local
    fi

    mkdir rd_plots/

    tar zxvf ~{rd_metrics_folder}
    Rscript ~{Rdtest_V2_script} \
      -b ~{rdtest_bed} \
      -m ~{rd_median_file} \
      -c rd_metrics_folder/ \
      -p TRUE \
      -o rd_plots/ \
      -s 10000000000 \
      ~{whitelist_arg}

    tar czvf ~{prefix}.tar.gz rd_plots

   >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File rd_plots = "~{prefix}.tar.gz"
  }

}

task RemoteTabixRdMetrics {
  input {
    File rdtest_bed
    File sample_batch_list
    String sv_pipeline_base_docker
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File rd_metrics_folder = "rd_metrics_folder.tar.gz"
    File remote_tabix_script = "remote_tabix.sh"
    File index_remote_tabix_script = "index_remote_tabix.sh"
  }
  command <<<

    set -euo pipefail

    cut -f5 ~{rdtest_bed} | sed -e 's/,/\n/g' | sort | uniq > sample_list
    grep -wf sample_list ~{sample_batch_list} | cut -f2 | sort | uniq > batch_list
    awk '{print $1, $2-100000, $3+100000}' ~{rdtest_bed} | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n | bedtools merge -i - > rdtest.bed
    mkdir rd_metrics_folder/
    
    #paste -d ' ' <(sed -e "s/^/    tabix -h -R /" batch_list | sed -e 's/$/ | bgzip -c >/') <(sed -e 's/\//\t/g' batch_list | awk '{print $NF}' | sed -e 's/^/rd_metrics_folder\//') | sed -e 's/^/GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` \ \\\n/' > remote_tabix.sh

    paste -d ' ' <(sed -e "s/^/tabix -h -R rdtest.bed /" batch_list | sed -e 's/$/ | bgzip -c >/') <(sed -e 's/\//\t/g' batch_list | awk '{print $NF}' | sed -e 's/^/rd_metrics_folder\//') | sed -e 's/^/GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` \ \\\n/' > remote_tabix.sh
    bash remote_tabix.sh

    ls rd_metrics_folder/*gz | sed -e 's/^/tabix -e 2 -b 2 /' > index_remote_tabix.sh
    bash index_remote_tabix.sh
    tar czvf rd_metrics_folder.tar.gz rd_metrics_folder/

 
   >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CollectRdMedian {
  input {
    File rdtest_bed
    File sample_batch_list
    String sv_pipeline_base_docker
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File rd_median_file = "rd_median.txt"
  }
  command <<<

    set -euo pipefail

    cut -f5 ~{rdtest_bed} | sed -e 's/,/\n/g' | sort | uniq > sample_list
    grep -wf sample_list ~{sample_batch_list} | cut -f3 | sort | uniq > batch_median_list
    mkdir rd_median_folder
    sed -e 's/^/gsutil cp /' batch_median_list | sed -e 's/$/ .\/rd_median_folder\//' > download_batch_median_list.sh
    bash download_batch_median_list.sh
    paste rd_median_folder/* > rd_median.txt
 
   >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


