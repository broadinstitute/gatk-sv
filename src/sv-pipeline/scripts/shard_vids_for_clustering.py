import argparse
import sys
import os
import pysam

def shard_vids_for_clustering(clustered_vcf, prefix, records_per_shard):
    vcf = pysam.VariantFile(clustered_vcf)

    # Exit early if the vcf is empty
    is_empty = True
    for record in vcf:
        is_empty = False
        break
    if is_empty:
        print("empty vcf - no shards will be produced")
        sys.exit(0)
    vcf.reset()

    current_cluster = None
    current_cluster_vids = []
    current_shard = 0
    current_shard_size = 0
    shard_path_format = f"{prefix}.vids.shard_{{}}.list"
    shard_path = shard_path_format.format(current_shard)
    fout = open(shard_path, 'w')
    if fout is None:
        raise IOError("Could not open '{}'".format(shard_path))
        # sys.exit(1)

    for record in vcf.fetch():
        cluster_id = record.info['CLUSTER']
        if cluster_id == current_cluster:
            current_cluster_vids.append(record.id)
        else:
            for vid in current_cluster_vids:
                fout.write(vid + '\n')
            current_shard_size += len(current_cluster_vids)
            if current_shard_size >= records_per_shard:
                current_shard += 1
                current_shard_size = 0
                fout.close()
                shard_path = shard_path_format.format(current_shard)
                fout = open(shard_path, 'w')
                if fout is None:
                    raise IOError("Could not open '{}'".format(shard_path))
                    # sys.exit(1)
            current_cluster_vids = [record.id]
            current_cluster = cluster_id

    # Write last cluster
    for vid in current_cluster_vids:
        fout.write(vid + '\n')
    current_shard_size += len(current_cluster_vids)
    fout.close()

    # Delete trailing empty shard
    if current_shard > 0 and current_shard_size == 0:
        os.remove(shard_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--clustered_vcf", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--records_per_shard", required=True)

    args = parser.parse_args()

    shard_vids_for_clustering(args.clustered_vcf, args.prefix, int(args.records_per_shard))
