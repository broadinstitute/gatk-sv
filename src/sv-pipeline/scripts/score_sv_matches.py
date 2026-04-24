import sys
import argparse
import pandas as pd
import pysam
import numpy as np


NON_INTERVAL_TYPES = {'CTX', 'BND'}
INS_CPX_TYPES = {'dDUP', 'dDUP_iDEL', 'INS_iDEL'}


def calculate_log_af_dif(af, truth_af):
    if af is None or truth_af is None:
        return None
    return abs(np.log10(af) - np.log10(truth_af))


def get_metric_from_list(r, metric, i):
    f = r.info.get(metric, None)
    if f is not None and isinstance(f, tuple):
        f = f[i]
    return f


def load_matches(vcf_path, this_cohort, that_cohort, cpx_intervals=False):
    print(f"Loading {vcf_path}")
    columns = ['VID_' + this_cohort, 'RECIPROCAL_OVERLAP', 'SIZE_SIMILARITY', 'BPDIST', 'logAF_DIF', 'VID_' + that_cohort, 'TYPE_PENALTY']
    dat = []
    with pysam.VariantFile(vcf_path) as vcf:
        for r in vcf:
            vid = r.id
            if cpx_intervals:
                vid = r.info.get('ORIGINAL_VID', r.id)
            svtype = r.info.get('SVTYPE')
            truth_vids = r.info.get('TRUTH_VID', None)
            cpxtype = r.info.get('CPX_TYPE', None)

            if truth_vids is not None:
                # only keep matches
                for i in range(len(truth_vids)):
                    # one line per variant pair
                    truth_svtype = get_metric_from_list(r, 'TRUTH_SVTYPE', i)
                    type_penalty = 0 if (svtype == truth_svtype) else 1  # penalize cross-type matches

                    # AF
                    # if comparing a multiallelic and a biallelic CNV, regardless of truth/eval, use RD_CN_ESTIMATED_AF
                    if svtype == "CNV" or truth_svtype == "CNV":
                        af = r.info.get('RD_CN_ESTIMATED_AF', None)
                    else:
                        af = r.info.get('AF', None)[0]

                    # TRUTH_AF
                    if svtype == "CNV" or truth_svtype == "CNV":
                        truth_af = get_metric_from_list(r, 'TRUTH_RD_CN_ESTIMATED_AF', i)
                    else:
                        truth_af = get_metric_from_list(r, 'TRUTH_AF', i)

                    # BKPT DIST
                    dist_start = get_metric_from_list(r, 'TRUTH_DISTANCE_START', i)
                    dist_end = get_metric_from_list(r, 'TRUTH_DISTANCE_END', i)
                    bd = max(dist_start, dist_end)
                    if svtype == 'CPX':
                        max_start = max(r.info.get('TRUTH_INTERVAL_START_DISTANCE', None))
                        max_end = max(r.info.get('TRUTH_INTERVAL_END_DISTANCE', None))
                        bd = max(bd, max_start, max_end)


                    # RO & SS
                    ro = get_metric_from_list(r, 'TRUTH_RECIPROCAL_OVERLAP', i)
                    ss = get_metric_from_list(r, 'TRUTH_SIZE_SIMILARITY', i)
                    if svtype in NON_INTERVAL_TYPES:
                        ro = 1
                        ss = 1
                    elif svtype == 'CPX':
                        interval_ro = r.info.get('TRUTH_INTERVAL_RECIPROCAL_OVERLAP', None)
                        interval_ss = r.info.get('TRUTH_INTERVAL_SIZE_SIMILARITY', None)
                        intervals = r.info.get('CPX_INTERVALS', None)
                        # for insertion-type CPX variants extract source & sink (but not cross-subtype matches)
                        if cpxtype in INS_CPX_TYPES and not cpx_intervals:
                            source_type = "DUP"
                            if cpxtype == "INS_iDEL":
                                source_type = "INS"
                            source_idx = [i for i in range(len(intervals)) if intervals[i].startswith(source_type)][0]
                            source_ro = interval_ro[source_idx]
                            source_ss = interval_ss[source_idx]
                            sink_ro = ro
                            sink_ss = ss
                            if cpxtype.endswith("iDEL"):
                                sink_idx = [i for i in range(len(intervals)) if intervals[i].startswith('DEL_')][0]
                                sink_ro = interval_ro[sink_idx]
                                sink_ss = interval_ss[sink_idx]
                            ro = (source_ro + sink_ro)/2
                            ss = (source_ss + sink_ss)/2
                        else:
                            # cross-subtype matches: always average across intervals - only relevant intervals included
                            ro = np.mean(interval_ro)
                            ss = np.mean(interval_ss)

                    # calculate and add BPDIST, logAF_DIF, TRUTH_VID as VID_A/B
                    r_data = [vid, ro, ss, bd, calculate_log_af_dif(af, truth_af), truth_vids[i], type_penalty]
                    dat.append(r_data)
    return pd.DataFrame(dat, columns=columns)



def merge_tables(df_a, df_b):
    # merge A-to-B and B-to-A tables
    merged = df_a.merge(df_b, on=['VID_A', 'VID_B'], how='outer')
    merged['RO'] = merged.RECIPROCAL_OVERLAP_x.where(merged.RECIPROCAL_OVERLAP_x.notnull(), merged.RECIPROCAL_OVERLAP_y)
    merged['SIZESIM'] = merged.SIZE_SIMILARITY_x.where(merged.SIZE_SIMILARITY_x.notnull(), merged.SIZE_SIMILARITY_y)
    merged['BPDIST'] = merged.BPDIST_x.where(merged.BPDIST_x.notnull(), merged.BPDIST_y)
    merged['logAF_DIF'] = merged.logAF_DIF_x.where(merged.logAF_DIF_x.notnull(), merged.logAF_DIF_y)

    # if logafdif is still null (RD_CN_ESTIMATED_AF can be null when no RD_CN), set to max
    max_log_af_dif = np.nanmax(merged.logAF_DIF.values)
    merged['logAF_DIF'] = merged.logAF_DIF.where(merged.logAF_DIF.notnull(), max_log_af_dif)

    # rank metrics
    merged['RANK_RO'] = merged.RO.rank()
    merged['RANK_AF_DIF'] = merged.logAF_DIF.rank(ascending=False)
    merged['RANK_BPDIST'] = merged.BPDIST.rank(ascending=False)
    merged['RANK_SS'] = merged.SIZESIM.rank()

    # calculate ranksum
    merged['ORIG_SCORE'] = merged.RANK_SS + merged.RANK_AF_DIF + merged.RANK_BPDIST

    # penalize cross-type matches by giving them negative scores, but retain order by score
    merged['SCORE'] = merged.ORIG_SCORE.where(merged.TYPE_PENALTY == 0, 0 - np.nanmax(merged.ORIG_SCORE) + merged.ORIG_SCORE)
    return merged[['VID_A', 'VID_B', 'SCORE']]


def add_cpx(merged, merged_cpx):
    print("Combining cross-subtype complex interval matches with all other matches")
    # force cross-subtype CPX matches to have lower scores than same-subtype CPX matches by making them negative
    # but retain ascending order by match quality
    merged_cpx['ORIG_SCORE'] = merged_cpx['SCORE']
    merged_cpx['SCORE'] = 0 - np.nanmax(merged_cpx.ORIG_SCORE) + merged_cpx.ORIG_SCORE

    return pd.concat([merged, merged_cpx[['VID_A', 'VID_B', 'SCORE']]])


def _parse_arguments(argv):
    parser = argparse.ArgumentParser(
        description="Create VID match and score table for input to SVFederate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--a-to-b", type=str, required=True,
                        help="MatchSVs VCF with cohort A as eval and cohort B as truth")
    parser.add_argument("--b-to-a", type=str, required=True,
                        help="MatchSVs VCF with cohort B as eval and cohort A as truth")
    parser.add_argument("--a-to-b-cpx", type=str, required=True,
                        help="MatchSVs complex intervals VCF with cohort A as eval and cohort B as truth")
    parser.add_argument("--b-to-a-cpx", type=str, required=True,
                        help="MatchSVs complex intervals VCF with cohort B as eval and cohort A as truth")
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

    print("Merging tables and calculating scores")
    merged = merge_tables(df_a, df_b)

    cpx_a = load_matches(args.a_to_b_cpx, "A", "B", cpx_intervals=True)
    cpx_b = load_matches(args.b_to_a_cpx, "B", "A", cpx_intervals=True)

    print("Merging complex intervals tables and calculating scores")
    merged_cpx = merge_tables(cpx_a, cpx_b)

    scores = add_cpx(merged, merged_cpx)

    print(f"Printing variant pairs and scores to {args.out}")
    scores.to_csv(args.out, sep='\t', header=True, index=False)


if __name__ == "__main__":
    main()
