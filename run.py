"""CuraCascade — Headless CLI & Pipeline Orchestration Interface.

Command-line entrypoint for querying clinical evidence inconsistencies,
verifying database integrity and cryptographic provenance, and launching
the interactive clinical audit dashboard.

Exit Status Codes:
    0 — Success / clean audit (no critical unpropagated Class A findings flagged)
    1 — Critical unpropagated Class A findings flagged in audit
    2 — Integrity or security verification failure
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Ensure CuraCascade root is on sys.path
_ROOT_DIR = Path(__file__).resolve().parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))

# Ensure robust UTF-8 handling on Windows consoles
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import config as C


# ==============================================================================
# Audit Field Definitions
# ==============================================================================

AUDIT_ALL_FIELDS: List[str] = [
    "finding_id",
    "audit_class",
    "severity",
    "gene_symbol",
    "disease_id",
    "disease_name",
    "clingen_classification",
    "clingen_date",
    "clingen_gcep",
    "clingen_report_url",
    "target_type",
    "target_id",
    "target_name",
    "target_status",
    "evidence_summary",
    "submitters_or_panel_info",
    "last_evaluated_date",
    "days_unresolved",
    "gnomad_af",
    "mcaf_threshold",
    "population_refuted_flag",
    "created_at",
]


# ==============================================================================
# Database Query Engine
# ==============================================================================

def query_audit_records(
    db_path: Path,
    genes: Optional[List[str]] = None,
    audit_class: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Query audit_findings table with parameterized filters and deterministic sorting.

    Args:
        db_path: Target SQLite database file.
        genes: Optional list of uppercase gene symbols.
        audit_class: Optional exact audit class ('Class A', 'Class B', 'Class C').

    Returns:
        List of matching finding dictionaries.
    """
    if not db_path.exists():
        raise FileNotFoundError(f"Database file not found: {db_path}")

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    try:
        cursor = conn.cursor()

        # Verify that audit_findings table exists
        cursor.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='audit_findings';"
        )
        if not cursor.fetchone():
            raise sqlite3.OperationalError(
                "Table 'audit_findings' does not exist in the database. Please seed data first."
            )

        fields_sql = ", ".join(AUDIT_ALL_FIELDS)
        query = f"SELECT {fields_sql} FROM audit_findings"
        conditions: List[str] = []
        params: List[Any] = []

        if audit_class:
            conditions.append("audit_class = ?")
            params.append(audit_class)

        if genes:
            placeholders = ", ".join("?" for _ in genes)
            conditions.append(f"UPPER(gene_symbol) IN ({placeholders})")
            params.extend(genes)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        query += """
        ORDER BY 
            CASE audit_class 
                WHEN 'Class A' THEN 1 
                WHEN 'Class B' THEN 2 
                WHEN 'Class C' THEN 3 
                ELSE 4 
            END,
            days_unresolved DESC,
            gene_symbol ASC,
            finding_id ASC;
        """

        cursor.execute(query, params)
        rows = [dict(r) for r in cursor.fetchall()]
        return rows
    finally:
        conn.close()


# ==============================================================================
# Presentation & Formatting Utilities
# ==============================================================================

def truncate_cell(val: Any, max_len: int) -> str:
    """Format and truncate a cell value to fit terminal column width."""
    if val is None:
        return ""
    text = str(val).strip()
    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


def format_table_report(
    findings: List[Dict[str, Any]],
    db_path: Path,
    gene_filter_str: str,
    class_filter_str: str,
) -> str:
    """Format findings into an enterprise-grade terminal report with aligned columns."""
    border = "=" * 165
    sub_border = "-" * 165

    class_a_count = sum(1 for r in findings if r.get("audit_class") == "Class A")
    class_b_count = sum(1 for r in findings if r.get("audit_class") == "Class B")
    class_c_count = sum(1 for r in findings if r.get("audit_class") == "Class C")
    distinct_genes = len(set(r.get("gene_symbol") for r in findings if r.get("gene_symbol")))

    lines: List[str] = [
        border,
        "🧬 CURACASCADE CLINICAL EVIDENCE INCONSISTENCY AUDIT REPORT",
        f"Target Database : {db_path.resolve()}",
        f"Gene Filter     : {gene_filter_str}",
        f"Class Filter    : {class_filter_str}",
        f"Summary Counts  : {len(findings):,} Total Discrepancies across {distinct_genes} Genes "
        f"(Class A: {class_a_count}, Class B: {class_b_count}, Class C: {class_c_count})",
        border,
    ]

    # Columns configuration: (key, header_title, width, alignment)
    columns: List[Tuple[str, str, int, str]] = [
        ("finding_id", "ID", 5, "right"),
        ("audit_class", "Class", 8, "left"),
        ("severity", "Severity", 9, "left"),
        ("gene_symbol", "Gene", 8, "left"),
        ("disease_name", "Disease Indication", 26, "left"),
        ("clingen_classification", "ClinGen", 9, "left"),
        ("target_type", "Target Layer", 16, "left"),
        ("target_name", "Target Asset / Variant", 34, "left"),
        ("target_status", "Target Status", 26, "left"),
        ("days_unresolved", "Lag (Days)", 10, "right"),
    ]

    header_parts: List[str] = []
    for _, title, width, align in columns:
        if align == "right":
            header_parts.append(f"{title:>{width}}")
        else:
            header_parts.append(f"{title:<{width}}")
    lines.append(" | ".join(header_parts))
    lines.append(sub_border)

    if not findings:
        lines.append("No inconsistency findings matched the specified criteria.")
    else:
        for r in findings:
            row_parts: List[str] = []
            for key, _, width, align in columns:
                val = r.get(key)
                cell_text = truncate_cell(val, width)
                if align == "right":
                    row_parts.append(f"{cell_text:>{width}}")
                else:
                    row_parts.append(f"{cell_text:<{width}}")
            lines.append(" | ".join(row_parts))

    lines.append(sub_border)
    if class_a_count > 0:
        lines.append(
            f"🚨 CRITICAL AUDIT ALERT: Flagged {class_a_count:,} unpropagated Class A findings "
            f"across {distinct_genes} genes!"
        )
        lines.append(
            "   Clinical Risk: Pathogenic/Likely Pathogenic ClinVar assertions persist in genes "
            "refuted/disputed by ClinGen GCEPs."
        )
    else:
        lines.append(
            "✅ CLEAN AUDIT: 0 unpropagated Class A findings flagged. "
            "No high-confidence ClinVar contradictions detected."
        )
    lines.append(border)

    return "\n".join(lines)


def format_json_report(findings: List[Dict[str, Any]]) -> str:
    """Format findings as indented JSON string."""
    return json.dumps(findings, indent=2, default=str)


def format_csv_report(findings: List[Dict[str, Any]]) -> str:
    """Format findings as standard CSV text."""
    output = io.StringIO()
    fields = AUDIT_ALL_FIELDS
    if findings:
        # Preserve all keys present in dict
        keys = list(findings[0].keys())
        writer = csv.DictWriter(output, fieldnames=keys)
        writer.writeheader()
        writer.writerows(findings)
    else:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
    return output.getvalue()


# ==============================================================================
# Subcommand Handlers
# ==============================================================================

def handle_audit(args: argparse.Namespace) -> int:
    """Execute the headless clinical audit query and format output.

    Returns:
        0 — If no Class A findings are flagged (success / clean audit)
        1 — If critical unpropagated Class A findings are flagged
        2 — If database or query execution error occurs
    """
    # 1. Determine target database
    target_db = Path(args.db) if getattr(args, "db", None) else C.DB_PATH

    # 2. Parse gene filters
    genes_list: Optional[List[str]] = None
    gene_display_str = "ALL"
    if args.genes:
        parsed_genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        if parsed_genes:
            genes_list = parsed_genes
            gene_display_str = ", ".join(genes_list)

    # 3. Parse inconsistency class filter
    raw_cls = (args.audit_class or "all").strip().lower()
    class_map = {
        "a": "Class A",
        "class a": "Class A",
        "b": "Class B",
        "class b": "Class B",
        "c": "Class C",
        "class c": "Class C",
        "all": None,
        "any": None,
        "*": None,
    }
    if raw_cls not in class_map:
        print(
            f"Error: Invalid class filter '{args.audit_class}'. Expected 'A', 'B', 'C', or 'all'.",
            file=sys.stderr,
        )
        return 2

    filter_class = class_map[raw_cls]
    class_display_str = filter_class if filter_class else "all"

    # 4. Query findings from database
    try:
        findings = query_audit_records(
            db_path=target_db,
            genes=genes_list,
            audit_class=filter_class,
        )
    except FileNotFoundError as fnf_err:
        print(f"Error: {fnf_err}", file=sys.stderr)
        return 2
    except sqlite3.OperationalError as op_err:
        print(f"Database Operational Error: {op_err}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"Unexpected Query Failure: {exc}", file=sys.stderr)
        return 2

    # 5. Format output
    fmt = (args.format or "table").strip().lower()
    if fmt == "table":
        output_content = format_table_report(
            findings=findings,
            db_path=target_db,
            gene_filter_str=gene_display_str,
            class_filter_str=class_display_str,
        )
    elif fmt == "json":
        output_content = format_json_report(findings)
    elif fmt == "csv":
        output_content = format_csv_report(findings)
    else:
        print(f"Error: Unsupported output format '{args.format}'.", file=sys.stderr)
        return 2

    # 6. Check for Class A findings
    class_a_findings = [r for r in findings if r.get("audit_class") == "Class A"]
    class_a_count = len(class_a_findings)

    # 7. Write to output file or stdout
    if args.output:
        out_path = Path(args.output)
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(output_content, encoding="utf-8")
            print(f"✅ Successfully wrote {len(findings):,} audit records to: {out_path.resolve()}")
            if class_a_count > 0:
                print(
                    f"🚨 [ALERT] {class_a_count:,} critical unpropagated Class A findings flagged.",
                    file=sys.stderr,
                )
            else:
                print("✅ [STATUS] Clean audit: 0 Class A findings flagged.")
        except Exception as write_err:
            print(f"Error writing to output file {out_path}: {write_err}", file=sys.stderr)
            return 2
    else:
        # Standard stdout emission
        if fmt == "table":
            print(output_content)
        elif fmt == "json":
            sys.stdout.write(output_content + "\n")
            if class_a_count > 0:
                sys.stderr.write(
                    f"[ALERT] {class_a_count:,} critical unpropagated Class A findings flagged.\n"
                )
        elif fmt == "csv":
            sys.stdout.write(output_content)
            if class_a_count > 0:
                sys.stderr.write(
                    f"[ALERT] {class_a_count:,} critical unpropagated Class A findings flagged.\n"
                )

    # 8. Return standard exit code:
    # 1 for critical unpropagated Class A findings flagged, 0 for success/clean audit
    if class_a_count > 0:
        return 1
    return 0


def handle_verify(args: argparse.Namespace) -> int:
    """Execute complete database integrity, domain range, and security verification.

    Returns:
        0 — If all integrity and provenance checks pass
        2 — If any integrity or provenance check fails
    """
    try:
        from pipeline.verify_integrity import (
            print_integrity_report,
            run_database_integrity_verification,
        )
    except ImportError as imp_err:
        print(f"Error importing verification suite: {imp_err}", file=sys.stderr)
        return 2

    databases_to_scan: List[Path] = []
    if getattr(args, "db", None):
        databases_to_scan.append(Path(args.db))
    else:
        showcase_db = C.DB_PATH
        full_db = _ROOT_DIR / "data" / "curacascade_full.db"
        databases_to_scan.append(showcase_db)
        if full_db.exists():
            databases_to_scan.append(full_db)

    overall_passed = True
    for i, db in enumerate(databases_to_scan):
        verify_provenance = (i == 0)
        try:
            report = run_database_integrity_verification(
                db_path=db,
                verify_provenance=verify_provenance,
            )
            print_integrity_report(report)
            if not report.all_passed:
                overall_passed = False
        except Exception as scan_err:
            print(f"Verification Error on {db}: {scan_err}", file=sys.stderr)
            overall_passed = False

    if overall_passed:
        return 0
    return 2


def handle_serve(args: argparse.Namespace) -> int:
    """Launch the interactive Streamlit dashboard on port 8501.

    Returns:
        0 — On clean exit or user cancellation
        1 — On failure to launch
    """
    app_file = _ROOT_DIR / "app.py"
    if not app_file.is_file():
        print(f"Error: app.py not found at {app_file}", file=sys.stderr)
        return 1

    port = getattr(args, "port", 8501) or 8501
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_file),
        "--server.port",
        str(port),
    ]

    if getattr(args, "host", None):
        cmd.extend(["--server.address", str(args.host)])

    if getattr(args, "headless", False):
        cmd.extend(["--server.headless", "true"])

    print("=" * 80)
    print("🚀 LAUNCHING CURACASCADE CLINICAL AUDIT DASHBOARD")
    print(f"Target Script : {app_file.resolve()}")
    print(f"Listening Port: {port}")
    print(f"Command       : {' '.join(cmd)}")
    print("=" * 80)

    try:
        result = subprocess.run(cmd)
        return result.returncode
    except KeyboardInterrupt:
        print("\n🛑 Streamlit dashboard terminated by user.")
        return 0
    except Exception as exc:
        print(f"Failed to launch Streamlit dashboard: {exc}", file=sys.stderr)
        return 1


def handle_sync(args: argparse.Namespace) -> int:
    """Inspect or ingest upstream deltas.
    
    Returns:
        0 - Success
        2 - Error
    """
    try:
        from pipeline.delta_watcher import check_upstream_deltas, apply_deltas
    except ImportError as imp_err:
        print(f"Error importing delta_watcher: {imp_err}", file=sys.stderr)
        return 2

    target_db = args.db if getattr(args, "db", None) else str(C.DB_PATH)
    
    try:
        if getattr(args, "check", False):
            check_upstream_deltas(target_db)
        elif getattr(args, "apply", False):
            apply_deltas(target_db)
        else:
            print("Please specify --check or --apply with the sync command.", file=sys.stderr)
            return 2
        return 0
    except Exception as exc:
        print(f"Sync error: {exc}", file=sys.stderr)
        return 2


# ==============================================================================
# CLI Argument Parser Construction
# ==============================================================================

def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argparse CLI with subcommands."""
    parser = argparse.ArgumentParser(
        prog="python run.py",
        description=(
            "🧬 CuraCascade — Clinical Evidence Downgrade Auditor & Pipeline CLI.\n"
            "Audit downstream evidence cascades, verify database security & integrity, "
            "or launch the interactive clinical dashboard."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        title="Subcommands",
        dest="subcommand",
        help="Subcommand to execute",
    )

    # -------------------------------------------------------------------------
    # 1. 'audit' Subcommand
    # -------------------------------------------------------------------------
    audit_parser = subparsers.add_parser(
        "audit",
        help="Audit cross-database discrepancies (ClinVar & PanelApp in downgraded genes).",
        description="Query and audit downstream evidence cascades with gene and class filters.",
    )
    audit_parser.add_argument(
        "--genes",
        type=str,
        default=None,
        help="Comma-separated list of gene symbols to filter (e.g. 'ACTN2,BARD1,SCN1B').",
    )
    audit_parser.add_argument(
        "--class",
        dest="audit_class",
        type=str,
        default="all",
        choices=["A", "B", "C", "all", "a", "b", "c", "ALL"],
        metavar="{A,B,C,all}",
        help="Filter by inconsistency class: 'A', 'B', 'C', or 'all' (default: all).",
    )
    audit_parser.add_argument(
        "--format",
        type=str,
        default="table",
        choices=["table", "json", "csv", "TABLE", "JSON", "CSV"],
        metavar="{table,json,csv}",
        help="Output format: 'table', 'json', or 'csv' (default: table).",
    )
    audit_parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional filepath to write results (writes to stdout if omitted).",
    )
    audit_parser.add_argument(
        "--db",
        type=str,
        default=None,
        help=f"Custom SQLite database path to query (default: {C.DB_PATH}).",
    )

    # -------------------------------------------------------------------------
    # 2. 'verify' Subcommand
    # -------------------------------------------------------------------------
    verify_parser = subparsers.add_parser(
        "verify",
        help="Run the complete database integrity and security verification suite.",
        description=(
            "Run complete database integrity and security verification suite "
            "(pipeline/verify_integrity.py): B-Tree integrity, foreign keys, "
            "domain boundaries, and cryptographic SHA-256 provenance."
        ),
    )
    verify_parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Custom database file path to verify (defaults to scanning showcase & full DBs).",
    )

    # -------------------------------------------------------------------------
    # 3. 'serve' Subcommand
    # -------------------------------------------------------------------------
    serve_parser = subparsers.add_parser(
        "serve",
        help="Launch the Streamlit dashboard on port 8501 (streamlit run app.py).",
        description="Launch the interactive Streamlit clinical curator dashboard on port 8501.",
    )
    serve_parser.add_argument(
        "--port",
        type=int,
        default=8501,
        help="Port to bind Streamlit server to (default: 8501).",
    )
    serve_parser.add_argument(
        "--host",
        type=str,
        default=None,
        help="Optional host address to bind server to (e.g. 'localhost', '0.0.0.0').",
    )
    serve_parser.add_argument(
        "--headless",
        action="store_true",
        default=False,
        help="Run Streamlit in headless server mode.",
    )

    # -------------------------------------------------------------------------
    # 4. 'sync' Subcommand
    # -------------------------------------------------------------------------
    sync_parser = subparsers.add_parser(
        "sync",
        help="Automated Upstream Delta Watcher & Continuous Sync Pipeline.",
        description="Inspect upstream API for pending reappraisal deltas or ingest incremental updates.",
    )
    sync_parser.add_argument(
        "--check",
        action="store_true",
        help="Dry run / metadata inspection of upstream ClinGen/ClinVar deltas.",
    )
    sync_parser.add_argument(
        "--apply",
        action="store_true",
        help="Ingest incremental deltas and update audit tables.",
    )
    sync_parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Custom SQLite database path to query.",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Main CLI entrypoint."""
    parser = build_parser()

    # If no arguments provided, display help
    if argv is None and len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 0

    args = parser.parse_args(argv)

    if not args.subcommand:
        parser.print_help(sys.stderr)
        return 0

    if args.subcommand == "audit":
        return handle_audit(args)
    elif args.subcommand == "verify":
        return handle_verify(args)
    elif args.subcommand == "serve":
        return handle_serve(args)
    elif args.subcommand == "sync":
        return handle_sync(args)
    else:
        parser.print_help(sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())
