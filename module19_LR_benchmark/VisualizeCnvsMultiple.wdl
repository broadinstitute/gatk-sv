version 1.0

import "VisualizeCnvs.wdl" as VisualizeCnvs
workflow VisualizeCnvsMultiple {
  input{
    # Note vcf will be faster
    Array[File] vcf_or_bed_list  # bed columns: chrom,start,end,name,svtype,samples

    Array[String] prefix_list
    
    File ped_file
    
    Int min_size
    String flags

    Array[File] median_files
    Array[File] rd_files

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_rdtest
  }

  scatter(i in range(length(prefix_list))){
  	call VisualizeCnvs.VisualizeCnvs as VisualizeCnvs{
  		input:
  			vcf_or_bed = vcf_or_bed_list[i],
  			prefix = prefix_list[i],
  			median_files = median_files,
  			rd_files = rd_files,
  			ped_file = ped_file,
  			min_size = min_size,
  			flags = flags,
  			sv_pipeline_docker = sv_pipeline_docker,
  			runtime_attr_rdtest = runtime_attr_rdtest

  	}
  }

  ouput{
  	Array[File] rdtest_plot_tars = VisualizeCnvs.rdtest_plots
  }
 }








  }

