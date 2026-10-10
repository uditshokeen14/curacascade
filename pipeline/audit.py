"""CuraCascade — Cross-Database Inconsistency Audit Engine.

Executes relational SQL joins and classification logic across ClinGen, ClinVar,
and Genomics England PanelApp to identify downstream evidence cascades:
- Class A: Pathogenic/Likely Pathogenic ClinVar assertions in Disputed or Refuted genes.
- Class B: Pathogenic/Likely Pathogenic ClinVar assertions in Limited-evidence genes.
- Class C: Green (diagnostic-grade) PanelApp gene entries in Disputed or Refuted genes.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd

import config as C
from database.connection import (
    ensure_database_ready,
    execute_many,
    execute_write,
    query_df,
)


def run_audit(
    db_path: Optional[Union[str, Path]] = None,
    export_csv: bool = True,
) -> Dict[str, Any]:
    """Execute the full cross-layer inconsistency audit and populate audit_findings.

    Args:
        db_path: Target SQLite database file path.
        export_csv: If True, writes results to data/processed/audit_findings.csv.

    Returns:
        Dictionary containing summary metrics and audit DataFrames.
    """
    ensure_database_ready(db_path=db_path)

    # Clean out previous audit results to ensure idempotency
    execute_write("DELETE FROM audit_findings;", db_path=db_path)

    # -------------------------------------------------------------------------
    # 1. SQL Query for Class A & Class B (ClinVar P/LP inconsistencies)
    # -------------------------------------------------------------------------
    # Class A: Disputed or Refuted genes + ClinVar P/LP
    # Class B: Limited evidence genes + ClinVar P/LP
    clinvar_audit_sql = """
    SELECT 
        cc.classification AS clingen_classification,
        cc.classification_date AS clingen_date,
        cc.gcep AS clingen_gcep,
        cc.report_url AS clingen_report_url,
        g.gene_symbol,
        d.disease_id,
        d.disease_name,
        cv.variation_id,
        cv.variant_name,
        cv.clinical_significance,
        cv.stars,
        cv.review_status,
        cv.num_submitters,
        cv.submitter_names,
        cv.last_evaluated,
        cv.phenotype_list
    FROM clingen_curations cc
    JOIN genes g ON cc.gene_symbol = g.gene_symbol
    JOIN diseases d ON cc.disease_id = d.disease_id
    JOIN clinvar_variants cv ON g.gene_symbol = cv.gene_symbol
    WHERE cv.is_plp = 1
      AND cc.classification IN ('Disputed', 'Refuted', 'Limited')
    ORDER BY 
        CASE cc.classification 
            WHEN 'Refuted' THEN 1 
            WHEN 'Disputed' THEN 2 
            WHEN 'Limited' THEN 3 
            ELSE 4 
        END,
        cv.stars DESC,
        g.gene_symbol ASC;
    """
    clinvar_rows = query_df(clinvar_audit_sql, db_path=db_path)

    findings_to_insert: List[Tuple[Any, ...]] = []
    today = datetime.now().date()

    for _, row in clinvar_rows.iterrows():
        classification = str(row["clingen_classification"]).strip()
        clingen_date_str = str(row["clingen_date"]) if pd.notna(row["clingen_date"]) else ""

        days_unresolved = 0
        if clingen_date_str:
            try:
                c_date = datetime.strptime(clingen_date_str.split()[0], "%Y-%m-%d").date()
                days_unresolved = max(0, (today - c_date).days)
            except ValueError:
                pass

        stars = int(row["stars"]) if pd.notna(row["stars"]) else 0
        star_str = f"{stars} star{'s' if stars != 1 else ''}"
        target_status = f"{row['clinical_significance']} ({star_str})"

        submitter_info = (
            f"{row['submitter_names']} ({row['num_submitters']} submitter{'s' if row['num_submitters'] != 1 else ''})"
            if pd.notna(row["submitter_names"])
            else f"{row['num_submitters']} submitter(s)"
        )

        if classification in ("Disputed", "Refuted"):
            audit_class = C.AUDIT_CLASS_A
            severity = "CRITICAL"
            evidence_summary = (
                f"ClinGen GCEP ({row['clingen_gcep'] or 'Consortium'}) classified the gene-disease "
                f"relationship as '{classification}' on {clingen_date_str}. ClinVar continues to report "
                f"this variant as {row['clinical_significance']} with {star_str} confidence."
            )
        else:  # Limited
            audit_class = C.AUDIT_CLASS_B
            severity = "HIGH"
            evidence_summary = (
                f"ClinGen GCEP ({row['clingen_gcep'] or 'Consortium'}) appraised gene evidence as 'Limited'. "
                f"Variants should not exceed VUS classification according to ACMG/ClinGen recommendations, "
                f"yet ClinVar asserts {row['clinical_significance']}."
            )

        findings_to_insert.append((
            audit_class,
            severity,
            row["gene_symbol"],
            row["disease_id"],
            row["disease_name"],
            classification,
            clingen_date_str,
            row["clingen_gcep"],
            row["clingen_report_url"],
            "ClinVar Variant",
            str(row["variation_id"]),
            str(row["variant_name"]),
            target_status,
            evidence_summary,
            submitter_info,
            str(row["last_evaluated"]) if pd.notna(row["last_evaluated"]) else "",
            days_unresolved,
        ))

    # -------------------------------------------------------------------------
    # 2. SQL Query for Class C (PanelApp Green entries in Disputed/Refuted genes)
    # -------------------------------------------------------------------------
    panelapp_audit_sql = """
    SELECT 
        cc.classification AS clingen_classification,
        cc.classification_date AS clingen_date,
        cc.gcep AS clingen_gcep,
        cc.report_url AS clingen_report_url,
        g.gene_symbol,
        d.disease_id,
        d.disease_name,
        pa.panel_id,
        pa.panel_name,
        pa.panel_version,
        pa.confidence,
        pa.mode_of_inheritance,
        pa.phenotypes
    FROM clingen_curations cc
    JOIN genes g ON cc.gene_symbol = g.gene_symbol
    JOIN diseases d ON cc.disease_id = d.disease_id
    JOIN panelapp_entries pa ON g.gene_symbol = pa.gene_symbol
    WHERE LOWER(pa.confidence) = 'green'
      AND cc.classification IN ('Disputed', 'Refuted')
    ORDER BY 
        CASE cc.classification 
            WHEN 'Refuted' THEN 1 
            WHEN 'Disputed' THEN 2 
            ELSE 3 
        END,
        pa.panel_name ASC;
    """
    panelapp_rows = query_df(panelapp_audit_sql, db_path=db_path)

    for _, row in panelapp_rows.iterrows():
        classification = str(row["clingen_classification"]).strip()
        clingen_date_str = str(row["clingen_date"]) if pd.notna(row["clingen_date"]) else ""

        days_unresolved = 0
        if clingen_date_str:
            try:
                c_date = datetime.strptime(clingen_date_str.split()[0], "%Y-%m-%d").date()
                days_unresolved = max(0, (today - c_date).days)
            except ValueError:
                pass

        target_status = f"Green (Diagnostic-Grade, v{row['panel_version'] or '1.0'})"
        panel_info = f"PanelApp ID {row['panel_id']} ({row['mode_of_inheritance'] or 'Unknown MOI'})"

        evidence_summary = (
            f"Genomics England PanelApp lists this gene as diagnostic-grade 'Green' for diagnostic reporting, "
            f"directly conflicting with ClinGen's '{classification}' validity determination made on {clingen_date_str}."
        )

        findings_to_insert.append((
            C.AUDIT_CLASS_C,
            "CRITICAL",
            row["gene_symbol"],
            row["disease_id"],
            row["disease_name"],
            classification,
            clingen_date_str,
            row["clingen_gcep"],
            row["clingen_report_url"],
            "PanelApp Panel",
            str(row["panel_id"]),
            str(row["panel_name"]),
            target_status,
            evidence_summary,
            panel_info,
            "",  # no variant evaluation date for panels
            days_unresolved,
        ))

    # -------------------------------------------------------------------------
    # 3. Insert findings into audit_findings table
    # -------------------------------------------------------------------------
    insert_sql = """
    INSERT INTO audit_findings (
        audit_class, severity, gene_symbol, disease_id, disease_name,
        clingen_classification, clingen_date, clingen_gcep, clingen_report_url,
        target_type, target_id, target_name, target_status, evidence_summary,
        submitters_or_panel_info, last_evaluated_date, days_unresolved
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    if findings_to_insert:
        execute_many(insert_sql, findings_to_insert, db_path=db_path)

    # -------------------------------------------------------------------------
    # 4. Compute Summary Metrics
    # -------------------------------------------------------------------------
    metrics = get_audit_metrics(db_path=db_path)
    findings_df = get_all_findings(db_path=db_path)

    if export_csv and not findings_df.empty:
        csv_path = C.PROCESSED_DATA_DIR / "audit_findings.csv"
        findings_df.to_csv(csv_path, index=False)

    return {
        "metrics": metrics,
        "findings_df": findings_df,
    }


def get_audit_metrics(db_path: Optional[Union[str, Path]] = None) -> Dict[str, int]:
    """Retrieve key numerical indicators of corpus integrity and audit findings.

    Args:
        db_path: Target database path.

    Returns:
        Dictionary of aggregated metrics.
    """
    counts_sql = """
    SELECT
        (SELECT COUNT(*) FROM genes) AS total_genes,
        (SELECT COUNT(*) FROM genes WHERE is_downgraded = 1) AS downgraded_genes,
        (SELECT COUNT(*) FROM clingen_curations) AS total_curations,
        (SELECT COUNT(*) FROM clinvar_variants) AS total_variants,
        (SELECT COUNT(*) FROM clinvar_variants WHERE is_plp = 1) AS total_plp_variants,
        (SELECT COUNT(*) FROM panelapp_entries) AS total_panel_entries,
        (SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class A') AS class_a_count,
        (SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class B') AS class_b_count,
        (SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class C') AS class_c_count,
        (SELECT COUNT(DISTINCT gene_symbol) FROM audit_findings WHERE audit_class = 'Class A') AS class_a_genes,
        (SELECT COUNT(DISTINCT gene_symbol) FROM audit_findings WHERE audit_class = 'Class B') AS class_b_genes,
        (SELECT COUNT(DISTINCT gene_symbol) FROM audit_findings WHERE audit_class = 'Class C') AS class_c_genes,
        (SELECT COUNT(*) FROM audit_findings) AS total_findings
    """
    df = query_df(counts_sql, db_path=db_path)
    if df.empty:
        return {}
    return {k: int(v) for k, v in df.iloc[0].to_dict().items()}


def _build_audit_where_clause(
    audit_class: Optional[str] = None,
    gene_symbol: Optional[str] = None,
    target_types: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
    prescriptive_action: Optional[str] = None,
) -> Tuple[str, List[Any]]:
    """Construct SQL WHERE clause and parameters for audit findings filtering."""
    conditions = []
    params: List[Any] = []

    if audit_class and audit_class != "All":
        conditions.append("audit_class = ?")
        params.append(audit_class)

    if gene_symbol and gene_symbol.strip():
        conditions.append("UPPER(gene_symbol) LIKE ?")
        params.append(f"%{gene_symbol.strip().upper()}%")

    if target_types:
        if isinstance(target_types, str) and target_types not in ("All", "All Layers"):
            conditions.append("target_type = ?")
            params.append(target_types)
        elif isinstance(target_types, (list, tuple)) and target_types:
            filtered_targets = [t for t in target_types if t not in ("All", "All Layers")]
            if filtered_targets:
                placeholders = ",".join("?" for _ in filtered_targets)
                conditions.append(f"target_type IN ({placeholders})")
                params.extend(filtered_targets)

    if prescriptive_action and prescriptive_action not in ("All", "All Actions"):
        act_str = str(prescriptive_action)
        if "Benign" in act_str or "BA1" in act_str:
            conditions.append("population_refuted_flag = 1")
        elif "Demote Panel" in act_str:
            conditions.append("(COALESCE(population_refuted_flag, 0) = 0) AND (target_type LIKE '%Panel%' OR audit_class = 'Class C') AND (clingen_classification IN ('Refuted', 'Disputed') OR audit_class = 'Class C')")
        elif "Amber Tier" in act_str:
            conditions.append("(COALESCE(population_refuted_flag, 0) = 0) AND (target_type LIKE '%Panel%' OR audit_class = 'Class C') AND (clingen_classification = 'Limited' OR audit_class = 'Class B')")
        elif "Retract P/LP" in act_str or "Reclassify to VUS" in act_str:
            conditions.append("(COALESCE(population_refuted_flag, 0) = 0) AND NOT (target_type LIKE '%Panel%' OR audit_class = 'Class C') AND (clingen_classification IN ('Refuted', 'Disputed') OR audit_class = 'Class A')")
        elif "Cap at VUS" in act_str:
            conditions.append("(COALESCE(population_refuted_flag, 0) = 0) AND NOT (target_type LIKE '%Panel%' OR audit_class = 'Class C') AND (clingen_classification = 'Limited' OR audit_class = 'Class B')")

    where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
    return where_clause, params


def _build_audit_order_clause(sort_by: Optional[str] = None) -> str:
    """Build SQL ORDER BY clause based on sort directive."""
    if not sort_by:
        return """
        ORDER BY 
            CASE audit_class 
                WHEN 'Class A' THEN 1 
                WHEN 'Class B' THEN 2 
                WHEN 'Class C' THEN 3 
                ELSE 4 
            END,
            days_unresolved DESC,
            gene_symbol ASC
        """
    s = str(sort_by).lower()
    if "longest" in s or ("latency" in s and "desc" in s) or "days ↓" in s:
        return "ORDER BY days_unresolved DESC, CASE audit_class WHEN 'Class A' THEN 1 WHEN 'Class B' THEN 2 WHEN 'Class C' THEN 3 ELSE 4 END, gene_symbol ASC"
    elif "shortest" in s or ("latency" in s and "asc" in s) or "days ↑" in s:
        return "ORDER BY days_unresolved ASC, CASE audit_class WHEN 'Class A' THEN 1 WHEN 'Class B' THEN 2 WHEN 'Class C' THEN 3 ELSE 4 END, gene_symbol ASC"
    elif "priority" in s or "class a" in s:
        return "ORDER BY CASE audit_class WHEN 'Class A' THEN 1 WHEN 'Class B' THEN 2 WHEN 'Class C' THEN 3 ELSE 4 END, days_unresolved DESC, gene_symbol ASC"
    elif "gene symbol (a" in s or ("gene" in s and "a → z" in s):
        return "ORDER BY gene_symbol ASC, days_unresolved DESC"
    elif "gene symbol (z" in s or ("gene" in s and "z → a" in s):
        return "ORDER BY gene_symbol DESC, days_unresolved DESC"
    elif "date" in s and ("newest" in s or "desc" in s):
        return "ORDER BY clingen_date DESC, days_unresolved DESC"
    elif "date" in s and ("oldest" in s or "asc" in s):
        return "ORDER BY clingen_date ASC, days_unresolved DESC"
    else:
        return """
        ORDER BY 
            CASE audit_class 
                WHEN 'Class A' THEN 1 
                WHEN 'Class B' THEN 2 
                WHEN 'Class C' THEN 3 
                ELSE 4 
            END,
            days_unresolved DESC,
            gene_symbol ASC
        """


def get_findings_count(
    db_path: Optional[Union[str, Path]] = None,
    audit_class: Optional[str] = None,
    gene_symbol: Optional[str] = None,
    target_types: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
    prescriptive_action: Optional[str] = None,
) -> int:
    """Retrieve total count of matching audit findings for pagination."""
    where_clause, params = _build_audit_where_clause(
        audit_class=audit_class,
        gene_symbol=gene_symbol,
        target_types=target_types,
        prescriptive_action=prescriptive_action,
    )
    sql = f"SELECT COUNT(*) FROM audit_findings {where_clause};"
    df = query_df(sql, params=params if params else None, db_path=db_path)
    if not df.empty:
        return int(df.iloc[0, 0])
    return 0


def get_all_findings(
    db_path: Optional[Union[str, Path]] = None,
    audit_class: Optional[str] = None,
    gene_symbol: Optional[str] = None,
    target_types: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
    prescriptive_action: Optional[str] = None,
    sort_by: Optional[str] = None,
    limit: Optional[int] = None,
    offset: Optional[int] = None,
) -> pd.DataFrame:
    """Retrieve audit findings with optional filtering, sorting, and server-side pagination.

    Args:
        db_path: Target database path.
        audit_class: Optional class filter ('Class A', 'Class B', 'Class C').
        gene_symbol: Optional gene symbol filter (case-insensitive substring match).
        target_types: Optional database layer filter ('ClinVar Variant', 'PanelApp Panel', or collection).
        prescriptive_action: Optional prescriptive action filter.
        sort_by: Optional sort order directive.
        limit: Optional maximum number of rows to retrieve (page size).
        offset: Optional row offset for server-side pagination.

    Returns:
        Pandas DataFrame of matching audit records.
    """
    where_clause, params = _build_audit_where_clause(
        audit_class=audit_class,
        gene_symbol=gene_symbol,
        target_types=target_types,
        prescriptive_action=prescriptive_action,
    )
    order_clause = _build_audit_order_clause(sort_by=sort_by)

    pagination_clause = ""
    if limit is not None and int(limit) > 0:
        pagination_clause = f"LIMIT {int(limit)}"
        if offset is not None and int(offset) > 0:
            pagination_clause += f" OFFSET {int(offset)}"

    sql = f"""
    SELECT 
        finding_id,
        audit_class,
        severity,
        gene_symbol,
        disease_name,
        clingen_classification,
        clingen_date,
        clingen_gcep,
        target_type,
        target_id,
        target_name,
        target_status,
        submitters_or_panel_info,
        last_evaluated_date,
        days_unresolved,
        evidence_summary,
        clingen_report_url,
        gnomad_af,
        mcaf_threshold,
        population_refuted_flag
    FROM audit_findings
    {where_clause}
    {order_clause}
    {pagination_clause};
    """
    return query_df(sql, params=params if params else None, db_path=db_path)


def get_gene_history(
    gene_symbol: str,
    db_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Retrieve full cross-database timeline and linked records for a specific gene.

    Args:
        gene_symbol: Gene symbol (e.g. 'MYLK2', 'MYH7').
        db_path: Target database path.

    Returns:
        Dictionary containing gene metadata, curations, variants, panel memberships, and findings.
    """
    symbol = gene_symbol.strip().upper()

    gene_df = query_df(
        "SELECT * FROM genes WHERE UPPER(gene_symbol) = ?",
        params=(symbol,),
        db_path=db_path,
    )
    if gene_df.empty:
        return {}

    curations_df = query_df(
        """
        SELECT cc.*, d.disease_name 
        FROM clingen_curations cc
        JOIN diseases d ON cc.disease_id = d.disease_id
        WHERE UPPER(cc.gene_symbol) = ?
        ORDER BY cc.classification_date DESC
        """,
        params=(symbol,),
        db_path=db_path,
    )

    variants_df = query_df(
        """
        SELECT cv.*, d.disease_name 
        FROM clinvar_variants cv
        LEFT JOIN diseases d ON cv.disease_id = d.disease_id
        WHERE UPPER(cv.gene_symbol) = ?
        ORDER BY cv.is_plp DESC, cv.stars DESC, cv.last_evaluated DESC
        """,
        params=(symbol,),
        db_path=db_path,
    )

    panels_df = query_df(
        """
        SELECT * FROM panelapp_entries 
        WHERE UPPER(gene_symbol) = ?
        ORDER BY 
            CASE LOWER(confidence) 
                WHEN 'green' THEN 1 
                WHEN 'amber' THEN 2 
                ELSE 3 
            END,
            panel_name ASC
        """,
        params=(symbol,),
        db_path=db_path,
    )

    findings_df = query_df(
        """
        SELECT * FROM audit_findings 
        WHERE UPPER(gene_symbol) = ?
        ORDER BY audit_class ASC
        """,
        params=(symbol,),
        db_path=db_path,
    )

    return {
        "gene": gene_df.iloc[0].to_dict(),
        "curations": curations_df,
        "variants": variants_df,
        "panels": panels_df,
        "findings": findings_df,
    }


def main():
    """CLI Entry point for running the audit pipeline."""
    parser = argparse.ArgumentParser(description="Run the CuraCascade cross-database audit.")
    parser.add_argument("--db", type=str, default=str(C.DB_PATH), help="Path to SQLite database.")
    args = parser.parse_args()

    print(f"Running CuraCascade Audit against: {args.db}")
    results = run_audit(db_path=args.db)
    metrics = results["metrics"]

    print("\n" + "=" * 60)
    print("  CURACASCADE AUDIT SUMMARY REPORT")
    print("=" * 60)
    print(f"  Total Genes Curated:           {metrics.get('total_genes', 0):>6,}")
    print(f"  Downgraded Genes (Exposure):   {metrics.get('downgraded_genes', 0):>6,}")
    print(f"  Total ClinVar Variants:        {metrics.get('total_variants', 0):>6,}")
    print(f"  Pathogenic/Likely Pathogenic:  {metrics.get('total_plp_variants', 0):>6,}")
    print(f"  Genomics England Panel Genes:  {metrics.get('total_panel_entries', 0):>6,}")
    print("-" * 60)
    print(f"  [CRITICAL] Class A Inconsistencies: {metrics.get('class_a_count', 0):>6,} ({metrics.get('class_a_genes', 0)} genes)")
    print(f"             (P/LP in Disputed/Refuted)")
    print(f"  [HIGH]     Class B Inconsistencies: {metrics.get('class_b_count', 0):>6,} ({metrics.get('class_b_genes', 0)} genes)")
    print(f"             (P/LP in Limited-Evidence)")
    print(f"  [CRITICAL] Class C Inconsistencies: {metrics.get('class_c_count', 0):>6,} ({metrics.get('class_c_genes', 0)} genes)")
    print(f"             (Diagnostic Green in Disputed/Refuted)")
    print("=" * 60)
    print(f"  Total Audit Findings Flagged:       {metrics.get('total_findings', 0):>6,}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
