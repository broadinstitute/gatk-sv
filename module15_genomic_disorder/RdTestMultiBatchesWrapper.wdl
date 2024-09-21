#script to generate rdplots across multiple batches
version 1.0

import "Structs.wdl"
import "RdTestMultiBatches.wdl" as rdtest_multi_batches
workflow RdTestMultiBatchesWrapper {
	input{
		Array[File] rdtest_bed_list
		Array[String] prefix_list
		File sample_batch_list
		File Rdtest_V2_script
		String sv_pipeline_base_docker
	}

	scatter (i in range(length(rdtest_bed_list))){
		call rdtest_multi_batches.RdTestMultiBatches{
			input:
				rdtest_bed = rdtest_bed_list[i],
				prefix = prefix_list[i],
				sample_batch_list = sample_batch_list,
				Rdtest_V2_script = Rdtest_V2_script,
				sv_pipeline_base_docker = sv_pipeline_base_docker
		}
	}

	output{
		Array[File] RdTestMultiBatchesList = RdTestMultiBatches.rd_plot_tarball
	}
	}