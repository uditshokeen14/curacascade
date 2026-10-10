import time
import pandas as pd
import sqlite3
import duckdb
from pathlib import Path

def benchmark_duckdb_vs_sqlite(db_path: str = None) -> dict:
    """Benchmark DuckDB columnar OLAP engine vs native SQLite."""
    if db_path:
        db_path = Path(db_path)
    else:
        db_path = Path("data/curacascade_full.db")
        if not db_path.exists():
            db_path = Path("data/curacascade.db")
            
    query = """
        SELECT 
            gnomad_popmax_ancestry,
            COUNT(*) as variant_count,
            AVG(gnomad_af_popmax / (mcaf_threshold + 1e-9)) as avg_fold_excess,
            SUM(population_refuted_flag) as refuted_count
        FROM clinvar_variants
        WHERE gnomad_popmax_ancestry IS NOT NULL
        GROUP BY gnomad_popmax_ancestry
        ORDER BY variant_count DESC;
    """
    
    # SQLite
    sqlite_start = time.perf_counter()
    try:
        sqlite_conn = sqlite3.connect(db_path)
        _ = pd.read_sql_query(query, sqlite_conn)
        sqlite_conn.close()
        sqlite_end = time.perf_counter()
        sqlite_time_ms = (sqlite_end - sqlite_start) * 1000
    except Exception as e:
        sqlite_time_ms = -1
        
    # DuckDB
    duckdb_start = time.perf_counter()
    try:
        duck_conn = duckdb.connect(":memory:")
        duck_conn.execute("INSTALL sqlite;")
        duck_conn.execute("LOAD sqlite;")
        duck_conn.execute(f"ATTACH '{db_path}' AS db (TYPE SQLITE);")
        # DuckDB query uses 'db.' prefix
        duck_query = query.replace("FROM clinvar_variants", "FROM db.clinvar_variants")
        _ = duck_conn.execute(duck_query).df()
        duck_conn.close()
        duckdb_end = time.perf_counter()
        duckdb_time_ms = (duckdb_end - duckdb_start) * 1000
    except Exception as e:
        duckdb_time_ms = -1

    speedup = sqlite_time_ms / duckdb_time_ms if duckdb_time_ms > 0 else 0
    
    return {
        "sqlite_ms": sqlite_time_ms,
        "duckdb_ms": duckdb_time_ms,
        "speedup_factor": speedup
    }

if __name__ == "__main__":
    results = benchmark_duckdb_vs_sqlite()
    print("OLAP Performance Benchmark")
    print("-" * 30)
    print(f"SQLite3 Time: {results['sqlite_ms']:.2f} ms")
    print(f"DuckDB Time:  {results['duckdb_ms']:.2f} ms")
    print(f"Speedup:      {results['speedup_factor']:.1f}x")
