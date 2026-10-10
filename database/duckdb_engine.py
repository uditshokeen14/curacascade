import duckdb
import pandas as pd
import sqlite3
import os

def get_duckdb_analytics(db_path: str = "data/curacascade.db") -> duckdb.DuckDBPyConnection:
    if not os.path.exists(db_path):
        db_path = "data/curacascade_full.db"
    
    conn = duckdb.connect(":memory:")
    # DuckDB's native SQLite scanner allows direct reading of SQLite files
    try:
        conn.execute("INSTALL sqlite;")
        conn.execute("LOAD sqlite;")
    except Exception:
        pass # May already be available
        
    conn.execute(f"ATTACH '{db_path}' AS db (TYPE SQLITE);")
    return conn

def get_fast_cohort_aggregates(db_path: str = "data/curacascade.db") -> pd.DataFrame:
    """Uses DuckDB to run multi-column GROUP BY queries returning Pandas DataFrames."""
    try:
        conn = get_duckdb_analytics(db_path)
        query = """
            SELECT 
                gnomad_popmax_ancestry,
                COUNT(*) as variant_count,
                AVG(gnomad_af_popmax / (mcaf_threshold + 1e-9)) as avg_fold_excess,
                SUM(population_refuted_flag) as refuted_count
            FROM db.clinvar_variants
            WHERE gnomad_popmax_ancestry IS NOT NULL
            GROUP BY gnomad_popmax_ancestry
            ORDER BY variant_count DESC;
        """
        df = conn.execute(query).df()
        conn.close()
        return df
    except Exception as e:
        print(f"DuckDB failed: {e}. Falling back to SQLite.")
        try:
            sqlite_conn = sqlite3.connect(db_path)
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
            df = pd.read_sql_query(query, sqlite_conn)
            sqlite_conn.close()
            return df
        except Exception as sqlite_e:
            print(f"SQLite fallback failed: {sqlite_e}")
            return pd.DataFrame()
