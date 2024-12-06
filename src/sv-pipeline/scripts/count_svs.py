import sys
from pathlib import Path

import duckdb


def validate_filters(con):
    sql = "SELECT svtype, min_svlen, max_svlen FROM sv_filters;"
    filters = con.sql(sql).fetchall()
    for f in filters:
        if f[1] < 0:
            raise ValueError("Min SV length must be >= 0")
        if f[1] > f[2]:
            raise ValueError("Min SV length must be <= max SV length")


def count_svs(con, filter_id):
    sql = """SELECT svtype, min_svlen, max_svlen
    FROM sv_filters
    WHERE id = ?;
    """
    filters = con.execute(sql, [filter_id]).fetchall()
    sql = (
        f"CREATE OR REPLACE TABLE sv_counts_{filter_id}"
        " AS SELECT sample, COUNT(*) AS count"
        " FROM svs_db.svs"
        " WHERE svtype = ? AND svlen >= ? AND svlen <= ?"
        " GROUP BY sample"
    )

    con.execute(sql, list(filters[0]))


counts_db = Path(sys.argv[1])
svs_db = Path(sys.argv[2])
if not counts_db.is_file():
    raise FileNotFoundError("Counts database not found")
if not svs_db.is_file():
    raise FileNotFoundError("SVs database not found")

with duckdb.connect(counts_db) as con:
    validate_filters(con)
    con.sql(f"ATTACH '{svs_db}' AS svs_db;")
    filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
    for i in filter_ids:
        count_svs(con, i[0])
