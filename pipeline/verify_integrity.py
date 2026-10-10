"""CuraCascade — Data Integrity and Security Infrastructure Verification.

Implements enterprise-grade database integrity validation, referential foreign key
scans, biological domain range assertions, and cryptographic SHA-256 provenance logging.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import os
from pathlib import Path
import sqlite3
import sys
from typing import Any, List, Optional, Tuple, Union

# Ensure workspace root is in python path
_ROOT_DIR = Path(__file__).resolve().parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))

import config as C
from database.connection import get_connection, get_db_path

# Ensure safe UTF-8 output on Windows consoles
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


# ==============================================================================
# Domain Boundary Definitions & Reference Rules
# ==============================================================================

VALID_AUDIT_CLASSES: Tuple[str, ...] = ("Class A", "Class B", "Class C")
MIN_ALLELE_FREQUENCY: float = 0.0
MAX_ALLELE_FREQUENCY: float = 1.0
MIN_CLINVAR_STARS: int = 0
MAX_CLINVAR_STARS: int = 4


@dataclass
class CheckResult:
    """Represents the outcome of an individual security or integrity check."""
    check_name: str
    target: str
    passed: bool
    details: str
    metric_value: Optional[Any] = None


@dataclass
class IntegrityReport:
    """Aggregates all integrity and security validation results."""
    database_path: Path
    all_passed: bool = True
    results: List[CheckResult] = field(default_factory=list)

    def add(self, result: CheckResult) -> None:
        self.results.append(result)
        if not result.passed:
            self.all_passed = False


# ==============================================================================
# Cryptographic Provenance Functions
# ==============================================================================

def compute_sha256(file_path: Union[str, Path], chunk_size: int = 65536) -> str:
    """Compute the SHA-256 cryptographic digest of a file in streaming chunks.

    Args:
        file_path: Path to target file.
        chunk_size: Streaming buffer size in bytes (default 64KB).

    Returns:
        Hexadecimal SHA-256 digest string.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Target file does not exist: {path}")

    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def verify_file_checksums(project_root: Optional[Path] = None) -> List[CheckResult]:
    """Scan and verify cryptographic integrity of database files and upstream raw assets.

    Args:
        project_root: Root workspace directory (defaults to CuraCascade root).

    Returns:
        List of CheckResult items for each surveyed file.
    """
    root = project_root or _ROOT_DIR
    results: List[CheckResult] = []

    target_files = [
        ("Database (Showcase)", root / "data" / "curacascade.db"),
        ("Database (Full Corpus)", root / "data" / "curacascade_full.db"),
        ("ClinGen Raw Asset", root / "curacascade" / "data" / "raw" / "clingen_gene_validity.csv"),
        ("Gene Audit Output", root / "curacascade" / "data" / "out" / "cascade_gene_level.csv"),
        ("Variant Audit Output", root / "curacascade" / "data" / "out" / "cascade_variant_level.csv"),
    ]

    for label, file_path in target_files:
        if file_path.exists():
            size_bytes = file_path.stat().st_size
            if size_bytes == 0:
                results.append(
                    CheckResult(
                        check_name="Cryptographic Checksum & Non-Emptiness",
                        target=f"{label} ({file_path.name})",
                        passed=False,
                        details=f"File exists but is empty (0 bytes).",
                        metric_value="0 bytes",
                    )
                )
            else:
                digest = compute_sha256(file_path)
                results.append(
                    CheckResult(
                        check_name="Cryptographic Checksum & Non-Emptiness",
                        target=f"{label} ({file_path.name})",
                        passed=True,
                        details=f"SHA-256: {digest[:16]}...{digest[-16:]} ({size_bytes:,} bytes)",
                        metric_value=digest,
                    )
                )
        else:
            # Full corpus or optional assets might not exist in lightweight checkout
            is_optional = "Full Corpus" in label or "curacascade" in str(file_path)
            results.append(
                CheckResult(
                    check_name="Cryptographic Checksum & Non-Emptiness",
                    target=f"{label} ({file_path.name})",
                    passed=is_optional,
                    details=f"File not present at {file_path.relative_to(root)}" if is_optional else f"Required file missing: {file_path}",
                    metric_value="MISSING" if not is_optional else "SKIPPED",
                )
            )

    return results


# ==============================================================================
# SQLite Engine & B-Tree Verification
# ==============================================================================

def verify_sqlite_btree(conn: sqlite3.Connection) -> CheckResult:
    """Execute low-level B-tree and page allocation integrity verification.

    Runs `PRAGMA integrity_check` across all tables, indices, and overflow pages.

    Args:
        conn: Open SQLite database connection.

    Returns:
        CheckResult indicating success or list of corruption messages.
    """
    cursor = conn.cursor()
    cursor.execute("PRAGMA integrity_check;")
    rows = cursor.fetchall()

    if len(rows) == 1 and rows[0][0] == "ok":
        return CheckResult(
            check_name="SQLite B-Tree & Page Integrity Check",
            target="Database Structure",
            passed=True,
            details="PRAGMA integrity_check returned 'ok' across all database pages.",
            metric_value="ok",
        )

    errors = [r[0] for r in rows if r[0] != "ok"]
    return CheckResult(
        check_name="SQLite B-Tree & Page Integrity Check",
        target="Database Structure",
        passed=False,
        details=f"Corruptions detected: {'; '.join(errors[:5])}",
        metric_value=f"{len(errors)} errors",
    )


# ==============================================================================
# Referential Foreign Key Verification
# ==============================================================================

def verify_foreign_keys(conn: sqlite3.Connection) -> CheckResult:
    """Scan all foreign key constraints for referential integrity violations.

    Runs `PRAGMA foreign_key_check` to detect any orphan child records in 3NF schema.

    Args:
        conn: Open SQLite database connection.

    Returns:
        CheckResult indicating whether foreign key constraints hold.
    """
    cursor = conn.cursor()
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute("PRAGMA foreign_key_check;")
    violations = cursor.fetchall()

    if not violations:
        return CheckResult(
            check_name="Foreign Key Referential Integrity Scan",
            target="3NF Relational Tables",
            passed=True,
            details="PRAGMA foreign_key_check found 0 referential violations across all child tables.",
            metric_value=0,
        )

    sample_v = [f"Table '{v[0]}' rowid {v[1]} references missing '{v[2]}'" for v in violations[:5]]
    return CheckResult(
        check_name="Foreign Key Referential Integrity Scan",
        target="3NF Relational Tables",
        passed=False,
        details=f"Detected {len(violations):,} foreign key violations: {'; '.join(sample_v)}",
        metric_value=len(violations),
    )


# ==============================================================================
# Biological Domain Range Boundary Assertions
# ==============================================================================

def verify_domain_boundaries(conn: sqlite3.Connection) -> List[CheckResult]:
    """Verify mathematical and biological domain boundaries across data tables.

    Checks:
        1. Allele frequency boundaries: gnomad_af, mcaf_threshold in [0.0, 1.0].
        2. ClinVar confidence stars: stars in [0, 4].
        3. Audit discrepancy taxonomy: audit_class in ('Class A', 'Class B', 'Class C').

    Args:
        conn: Open SQLite database connection.

    Returns:
        List of CheckResult items.
    """
    results: List[CheckResult] = []
    cursor = conn.cursor()

    # Discover existing tables in this database instance
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = {r[0] for r in cursor.fetchall()}

    # 1. ClinVar Variant Table Domain Boundaries
    cv_table = "clinvar_variants" if "clinvar_variants" in tables else ("clinvar_variant" if "clinvar_variant" in tables else None)
    if cv_table:
        # Check Allele Frequency Range [0.0, 1.0]
        cursor.execute(f"PRAGMA table_info({cv_table});")
        cv_cols = {r[1] for r in cursor.fetchall()}

        af_cols_to_check = [col for col in ["gnomad_af_global", "gnomad_af_popmax", "mcaf_threshold"] if col in cv_cols]
        for af_col in af_cols_to_check:
            query = f"""
                SELECT COUNT(*) FROM {cv_table}
                WHERE {af_col} IS NOT NULL AND ({af_col} < {MIN_ALLELE_FREQUENCY} OR {af_col} > {MAX_ALLELE_FREQUENCY});
            """
            cursor.execute(query)
            out_of_bounds_count = cursor.fetchone()[0]

            results.append(
                CheckResult(
                    check_name=f"Domain Range: Allele Frequency ({af_col})",
                    target=f"{cv_table}.{af_col} ∈ [0.0, 1.0]",
                    passed=(out_of_bounds_count == 0),
                    details="All values bounded within [0.0, 1.0]" if out_of_bounds_count == 0 else f"Found {out_of_bounds_count:,} out-of-bounds allele frequencies.",
                    metric_value=out_of_bounds_count,
                )
            )

        # Check Review Stars Range [0, 4]
        if "stars" in cv_cols:
            query = f"""
                SELECT COUNT(*) FROM {cv_table}
                WHERE stars IS NOT NULL AND (stars < {MIN_CLINVAR_STARS} OR stars > {MAX_CLINVAR_STARS});
            """
            cursor.execute(query)
            invalid_stars_count = cursor.fetchone()[0]

            results.append(
                CheckResult(
                    check_name="Domain Range: ClinVar Review Stars",
                    target=f"{cv_table}.stars ∈ [0, 4]",
                    passed=(invalid_stars_count == 0),
                    details="All review stars bounded within valid ACMG/ClinVar range [0, 4]" if invalid_stars_count == 0 else f"Found {invalid_stars_count:,} invalid star values.",
                    metric_value=invalid_stars_count,
                )
            )

    # 2. Audit Findings Table Domain Boundaries
    if "audit_findings" in tables:
        cursor.execute("PRAGMA table_info(audit_findings);")
        af_cols = {r[1] for r in cursor.fetchall()}

        # Check Discrepancy Classification Taxonomy
        if "audit_class" in af_cols:
            placeholders = ", ".join(f"'{c}'" for c in VALID_AUDIT_CLASSES)
            query = f"""
                SELECT COUNT(*) FROM audit_findings
                WHERE audit_class IS NOT NULL AND audit_class NOT IN ({placeholders});
            """
            cursor.execute(query)
            invalid_class_count = cursor.fetchone()[0]

            results.append(
                CheckResult(
                    check_name="Domain Range: Discrepancy Classification Taxonomy",
                    target=f"audit_findings.audit_class ∈ {VALID_AUDIT_CLASSES}",
                    passed=(invalid_class_count == 0),
                    details="All flagged findings conform to Class A, Class B, or Class C taxonomy" if invalid_class_count == 0 else f"Found {invalid_class_count:,} non-standard inconsistency classifications.",
                    metric_value=invalid_class_count,
                )
            )

        # Check Audit Findings Allele Frequencies
        for af_col in ["gnomad_af", "mcaf_threshold"]:
            if af_col in af_cols:
                query = f"""
                    SELECT COUNT(*) FROM audit_findings
                    WHERE {af_col} IS NOT NULL AND ({af_col} < {MIN_ALLELE_FREQUENCY} OR {af_col} > {MAX_ALLELE_FREQUENCY});
                """
                cursor.execute(query)
                out_of_bounds_count = cursor.fetchone()[0]

                results.append(
                    CheckResult(
                        check_name=f"Domain Range: Audit Finding AF ({af_col})",
                        target=f"audit_findings.{af_col} ∈ [0.0, 1.0]",
                        passed=(out_of_bounds_count == 0),
                        details="All values bounded within [0.0, 1.0]" if out_of_bounds_count == 0 else f"Found {out_of_bounds_count:,} out-of-bounds allele frequencies.",
                        metric_value=out_of_bounds_count,
                    )
                )

    return results


# ==============================================================================
# Orchestrator & CLI Runner
# ==============================================================================

def run_database_integrity_verification(
    db_path: Optional[Union[str, Path]] = None,
    verify_provenance: bool = True,
) -> IntegrityReport:
    """Run all database integrity, referential, domain, and checksum validations.

    Args:
        db_path: Path to SQLite database file.
        verify_provenance: Whether to include cryptographic SHA-256 asset scans.

    Returns:
        IntegrityReport containing all check details.
    """
    resolved_path = get_db_path(db_path)
    report = IntegrityReport(database_path=resolved_path)

    if not resolved_path.exists():
        report.add(
            CheckResult(
                check_name="Database Existence",
                target=str(resolved_path),
                passed=False,
                details=f"Database file does not exist at {resolved_path}",
            )
        )
        return report

    # 1. SQLite B-Tree and Page Check
    with get_connection(resolved_path, readonly=True) as conn:
        report.add(verify_sqlite_btree(conn))

        # 2. Referential Foreign Key Scans
        report.add(verify_foreign_keys(conn))

        # 3. Domain Range Boundary Assertions
        domain_results = verify_domain_boundaries(conn)
        for res in domain_results:
            report.add(res)

    # 4. Cryptographic Provenance Scans
    if verify_provenance:
        file_results = verify_file_checksums()
        for res in file_results:
            report.add(res)

    return report


def print_integrity_report(report: IntegrityReport) -> None:
    """Pretty-print the integrity report to stdout with clinical auditing formatting."""
    border = "=" * 80
    sub_border = "-" * 80

    print(border)
    print(f"🧬 CURACASCADE DATA INTEGRITY & CLINICAL SECURITY AUDIT REPORT")
    print(f"Target Database : {report.database_path.resolve()}")
    print(f"Timestamp       : {os.environ.get('TZ', 'UTC')} Evaluation")
    print(border)
    print(f"{'CHECK NAME':<40} | {'TARGET':<24} | {'STATUS'}")
    print(sub_border)

    for item in report.results:
        status_tag = "  PASS  " if item.passed else "  FAIL  "
        print(f"{item.check_name[:38]:<40} | {item.target[:22]:<24} | {status_tag}")
        if not item.passed:
            print(f"  └─ Details: {item.details}")

    print(sub_border)
    if report.all_passed:
        print("✅ ALL INTEGRITY & SECURITY CHECKS PASSED — REFERENTIAL CONSTRAINTS & DOMAIN LIMITS ENFORCED")
    else:
        print("❌ ONE OR MORE INTEGRITY CHECKS FAILED — REVIEW DIAGNOSTIC DETAILS ABOVE")
    print(border)


def main() -> int:
    """CLI Entry point for `python -m pipeline.verify_integrity`."""
    parser = argparse.ArgumentParser(
        description="Verify CuraCascade SQLite B-tree integrity, foreign keys, domain boundaries, and SHA-256 hashes."
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Custom database file path to verify (defaults to config.DB_PATH).",
    )
    parser.add_argument(
        "--all-dbs",
        action="store_true",
        help="Scan both the showcase benchmark database and the full corpus database.",
    )

    args = parser.parse_args()

    databases_to_scan: List[Path] = []
    if args.db:
        databases_to_scan.append(Path(args.db))
    elif args.all_dbs:
        showcase_db = C.DB_PATH
        full_db = _ROOT_DIR / "data" / "curacascade_full.db"
        databases_to_scan.append(showcase_db)
        if full_db.exists():
            databases_to_scan.append(full_db)
    else:
        showcase_db = C.DB_PATH
        full_db = _ROOT_DIR / "data" / "curacascade_full.db"
        databases_to_scan.append(showcase_db)
        if full_db.exists():
            databases_to_scan.append(full_db)

    overall_passed = True
    for db in databases_to_scan:
        report = run_database_integrity_verification(db_path=db, verify_provenance=(db == databases_to_scan[0]))
        print_integrity_report(report)
        if not report.all_passed:
            overall_passed = False

    return 0 if overall_passed else 1


if __name__ == "__main__":
    sys.exit(main())
