version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks


workflow PlotBenchmarkResults{
    input{
        File tp_query
        File tp_ref
        File fp_query
        File fp_ref 

        Boolean short_read_benchmark = false

        String sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.CalcuCompStat as calcu_comp_stat{
      input:
        tp_query = tp_query,
        tp_ref = tp_ref,
        fp_query = fp_query,
        fp_ref = fp_ref,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.PlotCompResults as plot_comp_results{
      input:
        comp_stat = calcu_comp_stat.comp_stat,
        docker_image = sv_pipeline_base_docker
    }

    output {
      File benchmark_stat = calcu_comp_stat.comp_stat
      File benchmark_figure = plot_comp_results.figure
  }
}

