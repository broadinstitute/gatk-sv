#!/bin/bash

import argparse
from pathlib import Path
import sys
from typing import List, Optional, Text


def write_script(
        bed_path: Text,
        reads_path: Text,
        snapshot_dir: Text,
        snapshot_prefix: Text,
        script_path: Text,
        reference_path: Text) -> None:
    Path(snapshot_dir).mkdir(parents=True, exist_ok=True)
    with open(bed_path) as bed, open(script_path, 'w') as script:
        script.write('new\n')
        script.write(f"genome {reference_path}\n")
        script.write(f"load {reads_path}\n")
        script.write(f"snapshotDirectory {snapshot_dir}\n")
        for line in bed:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            chrom = tokens[0]
            pos = int(tokens[1])
            end = int(tokens[2])
            vid = tokens[3]
            script.write(f"region {chrom} {pos} {end}\n")
            svsize = end - pos
            if end - pos < 10000:
                padding = max(svsize * 0.25, 500)
                script.write(f"goto {chrom}:{pos - padding}-{end + padding}\n")
                script.write('viewaspairs\n')
                script.write('collapse\n')
                script.write('sort insertsize\n')
                script.write(f"snapshot {snapshot_prefix}__{vid}.png\n")
            else:
                inner_padding = 500
                outer_padding = 1000
                # Start
                script.write(f"goto {chrom}:{pos - outer_padding}-{pos + inner_padding}\n")
                script.write('viewaspairs\n')
                script.write('collapse\n')
                script.write('sort insertsize\n')
                script.write(f"snapshot {snapshot_prefix}__{vid}__start.png\n")

                # End
                script.write(f"goto {chrom}:{end - inner_padding}-{end + outer_padding}\n")
                script.write('viewaspairs\n')
                script.write('collapse\n')
                script.write('sort insertsize\n')
                script.write(f"snapshot {snapshot_prefix}__{vid}__end.png\n")
        script.write('exit\n')


def parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create IGV screenshot script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bed", type=str, required=True,
                        help="Variant bed file (decompressed)")
    parser.add_argument("--reads", type=str, required=True,
                        help="BAM/CRAM URL or path")
    parser.add_argument("--snapshot-dir", type=str, required=True,
                        help="Output snapshot directory. It will be created if it doesn't exist.")
    parser.add_argument("--snapshot-prefix", type=str, required=True,
                        help="Output snapshot filename prefix")
    parser.add_argument("--out", type=str, required=True,
                        help="Output script path")
    parser.add_argument("--reference", type=str, required=True,
                        help="Reference fasta")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = parse_arguments(argv)
    write_script(
        bed_path=arguments.bed,
        reads_path=arguments.reads,
        snapshot_dir=arguments.snapshot_dir,
        snapshot_prefix=arguments.snapshot_prefix,
        script_path=arguments.out,
        reference_path=arguments.reference
    )


if __name__ == "__main__":
    main()
