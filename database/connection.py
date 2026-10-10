"""CuraCascade — Database Connection and Execution Manager.

Provides thread-safe connections, context managers, and high-performance
pandas/DuckDB bridge methods for relational bioinformatics data.
"""

from __future__ import annotations

import contextlib
import sqlite3
from pathlib import Path
from typing import Any, Generator, Iterable, Optional, Union

import pandas as pd

import config as C
from database.schema import init_db


def get_db_path(custom_path: Optional[Union[str, Path]] = None) -> Path:
    """Resolve the SQLite database file path.

    Args:
        custom_path: Optional custom path overriding default config.

    Returns:
        Resolved Path object.
    """
    if custom_path:
        return Path(custom_path)
    return C.DB_PATH


def ensure_genes_schema(conn: sqlite3.Connection) -> None:
    """Ensure that the genes table contains required schema columns."""
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('genes', 'gene');")
        tables = [r[0] for r in cursor.fetchall()]
        for tbl in tables:
            cursor.execute(f"PRAGMA table_info({tbl});")
            cols = [row[1] for row in cursor.fetchall()]
            if cols and "highest_downgrade_severity" not in cols:
                cursor.execute(f"ALTER TABLE {tbl} ADD COLUMN highest_downgrade_severity TEXT DEFAULT 'Disputed';")
                conn.commit()
    except Exception:
        pass


@contextlib.contextmanager
def get_connection(
    db_path: Optional[Union[str, Path]] = None,
    row_factory: bool = False,
    timeout: float = 30.0,
    readonly: bool = False,
) -> Generator[sqlite3.Connection, None, None]:
    """Yield a managed SQLite database connection.

    Enforces foreign keys, busy timeout, and WAL mode for high concurrency.
    When readonly=True, enforces PRAGMA query_only = ON to ensure frontend queries
    cannot inadvertently mutate relational database tables.

    Args:
        db_path: Target database path (defaults to config.DB_PATH).
        row_factory: If True, uses sqlite3.Row for dictionary-like column access.
        timeout: Lock acquisition timeout in seconds.
        readonly: If True, enforces PRAGMA query_only = ON or URI ?mode=ro.

    Yields:
        sqlite3.Connection instance.
    """
    path = get_db_path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not path.exists() and not readonly:
        pass
    elif not path.exists() and readonly:
        ensure_database_ready(path)

    if readonly:
        try:
            conn = sqlite3.connect(f"file:{path.resolve().as_posix()}?mode=ro", uri=True, timeout=timeout)
        except Exception:
            conn = sqlite3.connect(str(path), timeout=timeout)
            conn.execute("PRAGMA query_only = ON;")
    else:
        conn = sqlite3.connect(str(path), timeout=timeout)

    if row_factory:
        conn.row_factory = sqlite3.Row

    try:
        conn.execute("PRAGMA foreign_keys = ON;")
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        if readonly:
            conn.execute("PRAGMA query_only = ON;")
        else:
            ensure_genes_schema(conn)
        yield conn
    finally:
        conn.close()


def query_df(
    sql: str,
    params: Optional[Union[Iterable[Any], dict]] = None,
    db_path: Optional[Union[str, Path]] = None,
    readonly: bool = True,
) -> pd.DataFrame:
    """Execute a SELECT query and return results as a Pandas DataFrame.

    Hardened by default with PRAGMA query_only = ON and read-only URI connection mode.

    Args:
        sql: SQL query string.
        params: Optional query parameters (tuple or dict) for parameterized bindings.
        db_path: Optional database path override.
        readonly: Whether to enforce read-only connection (defaults to True).

    Returns:
        Pandas DataFrame of query results.
    """
    with get_connection(db_path, readonly=readonly) as conn:
        return pd.read_sql_query(sql, conn, params=params)


def execute_write(
    sql: str,
    params: Optional[Union[Iterable[Any], dict]] = None,
    db_path: Optional[Union[str, Path]] = None,
) -> int:
    """Execute a single INSERT/UPDATE/DELETE query with transaction commit.

    Args:
        sql: SQL statement string.
        params: Optional statement parameters for parameterized execution.
        db_path: Optional database path override.

    Returns:
        Number of rows affected.
    """
    with get_connection(db_path, readonly=False) as conn:
        cursor = conn.cursor()
        if params:
            cursor.execute(sql, params)
        else:
            cursor.execute(sql)
        conn.commit()
        return cursor.rowcount


def execute_many(
    sql: str,
    seq_of_params: Iterable[Union[Iterable[Any], dict]],
    db_path: Optional[Union[str, Path]] = None,
) -> int:
    """Execute batch INSERT/UPDATE statements within a transaction.

    Args:
        sql: SQL statement string.
        seq_of_params: Collection of parameter tuples or dicts.
        db_path: Optional database path override.

    Returns:
        Total number of executed parameter sets.
    """
    with get_connection(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(sql, seq_of_params)
        conn.commit()
        return cursor.rowcount


def ensure_database_ready(
    db_path: Optional[Union[str, Path]] = None,
    drop_existing: bool = False,
) -> None:
    """Ensure that the database file exists and has tables initialized.

    Args:
        db_path: Target database path.
        drop_existing: Whether to wipe and recreate tables.
    """
    with get_connection(db_path) as conn:
        init_db(conn, drop_existing=drop_existing)
