version 1.0

import "Structs.wdl"
import "CalculateGeneInterruptionPerGenomeUnit.wdl" as calculate_gene_interruption_per_genome_unit

workflow CalculateGeneInterruptionPerGenome {
    input{
        File src
        Array[File] AoU_bed_list
        Array[String] prefix_list
        Array[File] sample_list
        String sv_base_mini_docker
    }

    scatter(i in range(length(AoU_bed_list))){
            call calculate_gene_interruption_per_genome_unit.CalculateGeneInterruptionPerGenomeUnit as CalculateGeneInterruptionPerGenomeUnit{
            input:
                AoU_bed = AoU_bed_list[i],
                prefix = prefix_list[i],
                src = src,
                sample_list = sample_list,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    output{
        Array[File] SV_gene_sample_stat_list = CalculateGeneInterruptionPerGenomeUnit.SV_gene_sample_stat
    }
}


