import argparse
import os

# Implementation note:
# This code implements the same python code in the GetSampleBatchPEMap task without alterations.


def main(batch_sample_lists, batch_name_list, batch_pe_files, prefix):
    _batch_pe_files = []
    with open(batch_pe_files, 'r') as pe:
        for line in pe:
            local_file = os.path.basename(line.strip())
            _batch_pe_files.append(local_file)
    with open(f"{prefix}.sample_batch_pe_map.tsv", 'w') as out:
        for i in range(len(batch_name_list)):
            with open(batch_sample_lists[i], 'r') as inp:
                for line in inp:
                    sample = line.strip()
                    batch = batch_name_list[i]
                    pe_file = _batch_pe_files[i]
                    out.write(f"{sample}\t{batch}\t{pe_file}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-sample-lists", nargs='+', required=True)
    parser.add_argument("--batch-name-list", nargs='+', required=True)
    parser.add_argument("--batch-pe-files", required=True)
    parser.add_argument("--prefix", required=True)

    args = parser.parse_args()
    main(args.batch_sample_lists, args.batch_name_list, args.batch_pe_files, args.prefix)
