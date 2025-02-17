import argparse
from pathlib import Path

import duckdb


def dump_outliers(con, filter_id, outfile):
    sql = (
        "COPY (SELECT sample, count, svtype, min_svlen, max_svlen"
        " FROM ("
        f"SELECT sample, count, {filter_id} AS id"
        f" FROM outliers_{filter_id}"
        f" LEFT JOIN sv_counts_{filter_id}"
        " USING (sample))"
        " JOIN sv_filters USING (id))"
        f" TO '{outfile}'"
        " (DELIMITER '\t', HEADER false);"
    )
    con.sql(sql)


def dump(db, outdir):
    with duckdb.connect(db) as con:
        filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
        for i in filter_ids:
            outfile = Path(outdir, f"outliers_{i[0]}.tsv")
            dump_outliers(con, i[0], outfile)


def main(args):
    counts_db = Path(args.sv_counts_db)
    outdir = Path(args.outdir)
    if not counts_db.is_file():
        raise FileNotFoundError("Counts database not found")
    if outdir.is_dir():
        raise ValueError("Output directory exists")
    outdir.mkdir()

    dump(counts_db, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dump sample outliers from DuckDB file"
    )
    parser.add_argument(
        "sv_counts_db",
        metavar="SV_COUNTS_DB",
        help="Path to the input SV counts DuckDB datatbase",
    )
    parser.add_argument("outdir", metavar="OUTDIR", help="Path to the output directory")

    main(parser.parse_args())
