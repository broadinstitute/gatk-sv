version 1.0

import "AlignReads.wdl" as AR

workflow CallMinimap {
    input {
        File reads
        String map_preset
        String prefix
        File ref_map_file
        RuntimeAttr? runtime_attr_minimap2

    }

    parameter_meta {
        map_preset: "Options include: map-ont, map-hifi, map-pb, asm5, asm10, asm20"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call AR.Minimap2 as Align {
        input:
            reads = reads,
            ref_fasta = ref_map['fasta'],
            map_preset = map_preset,
            prefix = prefix,
            runtime_attr_override = runtime_attr_minimap2
    }

    output {
        File aligned_bam = Align.aligned_bam
        File aligned_bai = Align.aligned_bai
    }
}

