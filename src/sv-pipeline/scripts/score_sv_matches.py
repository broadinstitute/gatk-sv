import sys
import argparse
import pandas as pd
import pysam
import numpy as np


def calculate_log_af_dif(af, truth_af):
    if af is None or truth_af is None:
        return None
    return abs(np.log10(af) - np.log10(truth_af))


def load_matches(vcf_path, this_cohort, that_cohort):
    print(f"Loading {vcf_path}")
    # Info fields to load
    info_fields = ['TRUTH_RECIPROCAL_OVERLAP', 'TRUTH_SIZE_SIMILARITY', 'TRUTH_DISTANCE_START', 'TRUTH_DISTANCE_END', 'AF', 'TRUTH_AF']
    list_infos = {'AF'}
    columns = ['VID_' + this_cohort, 'TRUTH_RECIPROCAL_OVERLAP', 'TRUTH_SIZE_SIMILARITY', 'BPDIST', 'logAF_DIF', 'VID_' + that_cohort]
    dat = []
    with pysam.VariantFile(vcf_path) as vcf:
        for r in vcf:
            svtype = r.info.get('SVTYPE')
            filters = ",".join(f for f in r.filter)
            truth_vids = r.info.get('TRUTH_VID', None)
            if svtype in {'CPX', 'CTX'} or filters != "":
                continue  # ignore non-pass and variants without overlap metrics for now
            if truth_vids is not None:
                # only keep matches
                for i in range(len(truth_vids)):
                    # one line per variant pair
                    if "CNV" in truth_vids[i]:
                        continue  # ignore CNVs for now
                    r_data = [r.id]
                    dist_start = None
                    dist_end = None
                    af = None
                    truth_af = None
                    for k in info_fields:
                        f = r.info.get(k, None)
                        if k in list_infos:
                            f = f[0]
                        elif f is not None and isinstance(f, tuple):
                            f = f[i]
                        if k == 'TRUTH_DISTANCE_START':
                            dist_start = f
                            # store but don't add to the row data
                        elif k == 'AF':
                            af = f
                        elif k == 'TRUTH_AF':
                            truth_af = f
                        elif k == 'TRUTH_DISTANCE_END':
                            dist_end = f
                        else:
                            # record these values without modification
                            r_data.append(f)
                    # calculate and add BPDIST, logAF_DIF, TRUTH_VID as VID_A/B
                    r_data.extend([max(dist_start, dist_end), calculate_log_af_dif(af, truth_af), truth_vids[i]])
                    dat.append(r_data)
    return pd.DataFrame(dat, columns=columns)


def merge_tables(df_a, df_b):
    print("Merging tables and calculating scores")
    merged = df_a.merge(df_b, on=['VID_A', 'VID_B'], how='outer')
    merged['RO'] = merged.TRUTH_RECIPROCAL_OVERLAP_x.where(merged.TRUTH_RECIPROCAL_OVERLAP_x.notnull(), merged.TRUTH_RECIPROCAL_OVERLAP_y)
    # merged['SIZESIM'] = merged.TRUTH_SIZE_SIMILARITY_x.where(merged.TRUTH_SIZE_SIMILARITY_x.notnull(), merged.TRUTH_SIZE_SIMILARITY_y)
    merged['BPDIST'] = merged.BPDIST_x.where(merged.BPDIST_x.notnull(), merged.BPDIST_y)
    merged['logAF_DIF'] = merged.logAF_DIF_x.where(merged.logAF_DIF_x.notnull(), merged.logAF_DIF_y)
    merged['RANK_RO'] = merged.RO.rank()
    merged['RANK_AF_DIF'] = merged.logAF_DIF.rank(ascending=False)
    merged['RANK_BPDIST'] = merged.BPDIST.rank(ascending=False)
    # merged['RANK_SS'] = stats.rankdata(merged.SIZESIM)
    merged['SCORE'] = merged.RANK_RO + merged.RANK_AF_DIF + merged.RANK_BPDIST
    return merged[['VID_A', 'VID_B', 'SCORE']]


def _parse_arguments(argv):
    parser = argparse.ArgumentParser(
        description="Create VID match and score table for input to SVFederate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--a-to-b", type=str, required=True,
                        help="MatchSVs VCF with cohort A as eval and cohort B as truth")
    parser.add_argument("--b-to-a", type=str, required=True,
                        help="MatchSVs VCF with cohort B as eval and cohort A as truth")
    parser.add_argument("--out", type=str, required=True,
                        help="Output table")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main():
    args = _parse_arguments(sys.argv)

    df_a = load_matches(args.a_to_b, "A", "B")
    df_b = load_matches(args.b_to_a, "B", "A")

    merged = merge_tables(df_a, df_b)
    print(f"Printing variant pairs and scores to {args.out}")
    merged.to_csv(args.out, sep='\t', header=True, index=False)


if __name__ == "__main__":
    main()
