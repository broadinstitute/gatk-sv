#!/Users/kjaising/miniconda3/envs/venv/bin/python3

import argparse
import os
import json
import pandas as pd


def clean_vid(vid: str) -> str:
    if isinstance(vid, str) and vid.startswith('__') and vid.endswith('__'):
        return vid[2:-2]
    return vid


def load_lr_query_bed(file_path: str) -> pd.DataFrame:
    cols = ['chr', 'start', 'end', 'VID', 'svtype', 'length', 'AF', 'ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']
    return pd.read_csv(file_path, sep='\t', comment='#', header=None, names=cols)


def build_sample_result(sample_id: str, input_dirs: list, query_type: str = 'lr_query') -> dict:
    good_ids = set()
    bad_ids = set()

    for input_dir in input_dirs:
        query_file = os.path.join(input_dir, sample_id, f"{sample_id}.{query_type}.bed")
        if not os.path.exists(query_file):
            continue

        df = load_lr_query_bed(query_file)
        df['clean_VID'] = df['VID'].apply(clean_vid)

        good_mask = df['ovr1a'] != 'NO_OVR'
        good_ids.update(df.loc[good_mask, 'clean_VID'].astype(str).tolist())
        bad_ids.update(df.loc[~good_mask, 'clean_VID'].astype(str).tolist())

    bad_ids -= good_ids

    return {
        'good_variant_ids': sorted(good_ids),
        'bad_variant_ids': sorted(bad_ids)
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dirs', required=True, nargs='+', help='Directories with per-caller results (e.g., sniffles_results cutesv_results)')
    parser.add_argument('--mapping', required=True, help='TSV with SR_ID (and optional LR_ID) columns')
    parser.add_argument('--output-json', required=True, help='Path to write output JSON')
    parser.add_argument('--num-samples', type=int, default=None, help='Limit number of samples to include')
    parser.add_argument('--query-type', choices=['lr_query'], default='lr_query', help='Query type to read (fixed to lr_query)')

    args = parser.parse_args()

    mapping_df = pd.read_csv(args.mapping, sep='\t')
    samples = mapping_df['SR_ID'].tolist()

    results = {}
    processed = 0

    for sr_id in samples:
        if args.num_samples is not None and processed >= args.num_samples:
            break
        sample_result = build_sample_result(sr_id, args.input_dirs, args.query_type)
        results[sr_id] = sample_result
        processed += 1

    os.makedirs(os.path.dirname(args.output_json) or '.', exist_ok=True)
    with open(args.output_json, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Wrote JSON for {len(results)} samples to {args.output_json}")


if __name__ == '__main__':
    main() 