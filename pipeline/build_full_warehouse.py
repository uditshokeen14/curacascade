"""CuraCascade — Full-Scale Data Warehouse Builder (Objective 8).

Builds data/curacascade_full.db by ingesting full ClinGen gene-disease validity
dumps (~3,000+ curations), the full compiled ClinVar variant corpus (50,000+ variants),
and Genomics England PanelApp memberships (480+ entries).

Runs high-performance cross-layer relational joins to materialize the complete
inconsistency universe (Classes A, B, and C) with transactional batching and
B-tree indexing in < 2 minutes.
"""

from __future__ import annotations

import datetime
from pathlib import Path
import sqlite3
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8")
import time
from typing import Any, Dict, Optional, Union

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd
from pipeline.allele_frequency import (
    DISEASE_PREVALENCE_REGISTRY,
    GNOMAD_VARIANT_CATALOG,
    enrich_database_population_genetics,
)

INTERIM_DIR = _REPO_ROOT / "curacascade" / "data" / "interim"
OUT_DIR = _REPO_ROOT / "curacascade" / "data" / "out"
FULL_DB_PATH = _REPO_ROOT / "data" / "curacascade_full.db"
COMPAT_FULL_DB = OUT_DIR / "curacascade.db"

DOWNGRADED_CLASSIFICATIONS = {"Limited", "Disputed", "Refuted"}
SEVERITY_CLASS_A = {"Disputed", "Refuted"}


def build_full_warehouse(
    target_db_path: Optional[Union[str, Path]] = None,
    sync_compat: bool = True,
) -> Dict[str, Any]:
    """Execute the full-scale ETL pipeline and build data/curacascade_full.db."""
    start_time = time.time()
    db_path = Path(target_db_path or FULL_DB_PATH)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"🚀 [ETL] Starting CuraCascade Full Corpus Warehouse build -> {db_path.name}")

    # 1. Load Interim CSVs
    cg_path = INTERIM_DIR / "clingen_validity.csv.gz"
    cv_path = INTERIM_DIR / "clinvar_latest.csv.gz"
    pa_path = INTERIM_DIR / "panelapp_membership.csv.gz"

    if not cg_path.exists():
        raise FileNotFoundError(f"ClinGen validity dump not found: {cg_path}")

    print("   [1/6] Ingesting ClinGen gene-disease validity dump...")
    cg = pd.read_csv(cg_path)
    # Clean preamble / separators
    cg = cg[cg["gene"].notna() & (cg["gene"] != "+++++++++++")].copy()
    cg["gene"] = cg["gene"].astype(str).str.strip()
    cg["classification"] = cg["classification"].astype(str).str.strip()
    cg["classification_date"] = pd.to_datetime(cg["classification_date"], errors="coerce").dt.strftime("%Y-%m-%d")

    print(f"         Ingested {len(cg):,} ClinGen curations across {cg['gene'].nunique():,} unique genes.")

    print("   [2/6] Ingesting compiled ClinVar variant corpus...")
    cv = pd.read_csv(cv_path) if cv_path.exists() else pd.DataFrame()
    if not cv.empty:
        cv["gene"] = cv["gene"].astype(str).str.strip()
        cv["ClinicalSignificance"] = cv["ClinicalSignificance"].fillna("Uncertain significance")
        # Ensure stars column
        if "stars" not in cv.columns:
            cv["stars"] = cv["ReviewStatus"].apply(
                lambda s: 3 if "expert" in str(s).lower() else (2 if "multiple" in str(s).lower() else 1)
            )
        print(f"         Ingested {len(cv):,} ClinVar variants across {cv['gene'].nunique():,} genes.")

    print("   [3/6] Ingesting PanelApp panel memberships...")
    pa = pd.read_csv(pa_path) if pa_path.exists() else pd.DataFrame()
    if not pa.empty:
        pa["gene"] = pa["gene"].astype(str).str.strip()
        pa["confidence"] = pa["confidence"].astype(str).str.lower().str.strip()
        print(f"         Ingested {len(pa):,} PanelApp memberships across {pa['gene'].nunique():,} genes.")

    # 2. Database Initialization
    if db_path.exists():
        try:
            db_path.unlink()
        except Exception:
            pass

    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA synchronous = OFF;")
    conn.execute("PRAGMA journal_mode = MEMORY;")
    conn.execute("PRAGMA temp_store = MEMORY;")

    # 3. Create Schema Tables & Views
    conn.executescript(
        """
        DROP TABLE IF EXISTS audit_findings;
        DROP TABLE IF EXISTS variant_ancestry_metrics;
        DROP TABLE IF EXISTS downgrade_risk_predictions;
        DROP TABLE IF EXISTS panelapp_entries;
        DROP TABLE IF EXISTS clinvar_variants;
        DROP TABLE IF EXISTS clingen_curations;
        DROP TABLE IF EXISTS diseases;
        DROP TABLE IF EXISTS genes;

        CREATE TABLE genes (
            gene_symbol                TEXT PRIMARY KEY,
            hgnc_id                    TEXT,
            gene_name                  TEXT,
            chromosome                 TEXT DEFAULT 'N/A',
            is_downgraded              INTEGER NOT NULL DEFAULT 0,
            highest_downgrade_severity TEXT DEFAULT 'None'
        );

        CREATE TABLE diseases (
            disease_id             TEXT PRIMARY KEY,
            disease_name           TEXT NOT NULL,
            mondo_id               TEXT,
            prevalence_per_million REAL DEFAULT 1000.0,
            orphanet_id            TEXT
        );

        CREATE TABLE clingen_curations (
            curation_id         INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol         TEXT NOT NULL REFERENCES genes(gene_symbol),
            disease_id          TEXT REFERENCES diseases(disease_id),
            disease_name        TEXT,
            mondo_id            TEXT,
            moi                 TEXT,
            classification      TEXT NOT NULL,
            classification_date TEXT,
            gcep                TEXT,
            report_url          TEXT
        );

        CREATE TABLE clinvar_variants (
            variation_id          INTEGER PRIMARY KEY,
            gene_symbol           TEXT NOT NULL REFERENCES genes(gene_symbol),
            variant_name          TEXT,
            variant_type          TEXT,
            clinical_significance TEXT,
            review_status         TEXT,
            stars                 INTEGER DEFAULT 0,
            num_submitters        INTEGER DEFAULT 1,
            last_evaluated        TEXT,
            phenotype_list        TEXT,
            is_plp                INTEGER DEFAULT 1,
            disease_id            TEXT
        );

        CREATE TABLE panelapp_entries (
            entry_id            INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol         TEXT NOT NULL REFERENCES genes(gene_symbol),
            panel_id            INTEGER,
            panel_name          TEXT,
            panel_version       TEXT,
            confidence          TEXT,
            mode_of_inheritance TEXT,
            phenotypes          TEXT
        );

        CREATE TABLE audit_findings (
            finding_id               INTEGER PRIMARY KEY AUTOINCREMENT,
            audit_class              TEXT NOT NULL,
            severity                 TEXT NOT NULL,
            gene_symbol              TEXT NOT NULL,
            disease_id               TEXT,
            disease_name             TEXT,
            clingen_classification   TEXT NOT NULL,
            clingen_date             TEXT,
            clingen_gcep             TEXT,
            clingen_report_url       TEXT,
            target_type              TEXT NOT NULL,
            target_id                TEXT NOT NULL,
            target_name              TEXT NOT NULL,
            target_status            TEXT NOT NULL,
            evidence_summary         TEXT,
            submitters_or_panel_info TEXT,
            last_evaluated_date      TEXT,
            days_unresolved          INTEGER DEFAULT 0,
            gnomad_af                REAL DEFAULT 0.0,
            mcaf_threshold           REAL DEFAULT 0.0,
            population_refuted_flag  INTEGER DEFAULT 0,
            created_at               TEXT DEFAULT (datetime('now'))
        );
        """
    )

    # 4. Populate Core Entities
    print("   [4/6] Materializing relational entities (genes, diseases, curations)...")
    downgraded_gene_set = set(cg.loc[cg["classification"].isin(DOWNGRADED_CLASSIFICATIONS), "gene"])

    # Genes table
    genes_df = (
        cg.groupby("gene", as_index=False)
        .agg(hgnc_id=("hgnc_id", "first"))
        .rename(columns={"gene": "gene_symbol"})
    )
    genes_df["gene_name"] = genes_df["gene_symbol"]
    genes_df["chromosome"] = "N/A"
    genes_df["is_downgraded"] = genes_df["gene_symbol"].isin(downgraded_gene_set).astype(int)
    genes_df["highest_downgrade_severity"] = "None"
    limited_genes = set(cg.loc[cg["classification"] == "Limited", "gene"])
    disputed_genes = set(cg.loc[cg["classification"] == "Disputed", "gene"])
    refuted_genes = set(cg.loc[cg["classification"] == "Refuted", "gene"])
    genes_df.loc[genes_df["gene_symbol"].isin(limited_genes), "highest_downgrade_severity"] = "Limited"
    genes_df.loc[genes_df["gene_symbol"].isin(disputed_genes), "highest_downgrade_severity"] = "Disputed"
    genes_df.loc[genes_df["gene_symbol"].isin(refuted_genes), "highest_downgrade_severity"] = "Refuted"
    genes_df.to_sql("genes", conn, if_exists="append", index=False)

    # Diseases table
    disease_rows = []
    seen_diseases = set()
    for _, r in cg.iterrows():
        mid = str(r.get("mondo_id") or "").strip()
        dname = str(r.get("disease") or "").strip()
        did = mid if mid and mid != "nan" else dname
        if did and did not in seen_diseases:
            seen_diseases.add(did)
            prev = 1000.0
            orph = None
            if did in DISEASE_PREVALENCE_REGISTRY:
                prev = DISEASE_PREVALENCE_REGISTRY[did]["prevalence_per_million"]
                orph = DISEASE_PREVALENCE_REGISTRY[did]["orphanet_id"]
            disease_rows.append({
                "disease_id": did,
                "disease_name": dname or did,
                "mondo_id": mid if mid != "nan" else None,
                "prevalence_per_million": prev,
                "orphanet_id": orph,
            })
    if disease_rows:
        pd.DataFrame(disease_rows).to_sql("diseases", conn, if_exists="append", index=False)

    # ClinGen Curations
    cur_df = cg.rename(columns={"gene": "gene_symbol"}).copy()
    cur_df["disease_name"] = cur_df["disease"]
    cur_df["disease_id"] = cur_df["mondo_id"].fillna(cur_df["disease_name"])
    cols_cur = [
        "gene_symbol", "disease_id", "disease_name", "mondo_id",
        "moi", "classification", "classification_date", "gcep", "report_url"
    ]
    cur_df[[c for c in cols_cur if c in cur_df.columns]].to_sql(
        "clingen_curations", conn, if_exists="append", index=False
    )

    # ClinVar Variants
    if not cv.empty:
        cv_clean = cv.rename(
            columns={
                "gene": "gene_symbol",
                "VariationID": "variation_id",
                "Name": "variant_name",
                "Type": "variant_type",
                "ClinicalSignificance": "clinical_significance",
                "ReviewStatus": "review_status",
                "NumberSubmitters": "num_submitters",
                "LastEvaluated": "last_evaluated",
                "PhenotypeList": "phenotype_list",
            }
        ).drop_duplicates("variation_id")

        cv_clean["is_plp"] = cv_clean["clinical_significance"].apply(
            lambda s: 1 if "pathogenic" in str(s).lower() and "conflicting" not in str(s).lower() else 0
        )
        cv_clean["disease_id"] = None
        cols_cv = [
            "variation_id", "gene_symbol", "variant_name", "variant_type",
            "clinical_significance", "review_status", "stars", "num_submitters",
            "last_evaluated", "phenotype_list", "is_plp", "disease_id"
        ]
        cv_clean = cv_clean[cv_clean["gene_symbol"].isin(set(genes_df["gene_symbol"]))]
        cv_clean[[c for c in cols_cv if c in cv_clean.columns]].to_sql(
            "clinvar_variants", conn, if_exists="append", index=False
        )

    # PanelApp Entries
    if not pa.empty:
        pa_clean = pa.rename(columns={"gene": "gene_symbol"}).copy()
        cols_pa = [
            "gene_symbol", "panel_id", "panel_name", "panel_version",
            "confidence", "mode_of_inheritance", "phenotypes"
        ]
        pa_clean = pa_clean[pa_clean["gene_symbol"].isin(set(genes_df["gene_symbol"]))]
        pa_clean[[c for c in cols_pa if c in pa_clean.columns]].to_sql(
            "panelapp_entries", conn, if_exists="append", index=False
        )

    # Ingest Showcase Benchmark Variants & Panels for Complete Coverage
    showcase_db = Path("data/curacascade.db")
    if showcase_db.exists():
        with sqlite3.connect(str(showcase_db)) as s_conn:
            try:
                sc_cv = pd.read_sql_query("SELECT * FROM clinvar_variants", s_conn)
                if not sc_cv.empty:
                    existing_vids = set(pd.read_sql_query("SELECT variation_id FROM clinvar_variants", conn)["variation_id"])
                    sc_cv_new = sc_cv[~sc_cv["variation_id"].isin(existing_vids)].copy()
                    if not sc_cv_new.empty:
                        cols_cv_target = [
                            "variation_id", "gene_symbol", "variant_name", "variant_type",
                            "clinical_significance", "review_status", "stars", "num_submitters",
                            "last_evaluated", "phenotype_list", "is_plp", "disease_id"
                        ]
                        sc_cv_new[[c for c in cols_cv_target if c in sc_cv_new.columns]].to_sql(
                            "clinvar_variants", conn, if_exists="append", index=False
                        )
            except Exception as e:
                print(f"   [Note] Showcase ClinVar merge skipped: {e}")

            try:
                sc_pa = pd.read_sql_query("SELECT * FROM panelapp_entries", s_conn)
                if not sc_pa.empty:
                    cols_pa_target = [
                        "gene_symbol", "panel_id", "panel_name", "panel_version",
                        "confidence", "mode_of_inheritance", "phenotypes"
                    ]
                    existing_pa = set(tuple(x) for x in pd.read_sql_query("SELECT gene_symbol, panel_id FROM panelapp_entries", conn).values)
                    sc_pa_new = sc_pa[~sc_pa.apply(lambda r: (r["gene_symbol"], r["panel_id"]) in existing_pa, axis=1)].copy()
                    if not sc_pa_new.empty:
                        sc_pa_new[[c for c in cols_pa_target if c in sc_pa_new.columns]].to_sql(
                            "panelapp_entries", conn, if_exists="append", index=False
                        )
            except Exception as e:
                print(f"   [Note] Showcase PanelApp merge skipped: {e}")

    # 5. Cross-Layer Inconsistency Audit Materialization
    print("   [5/6] Executing cross-layer relational join to materialize Audit Findings (Classes A, B, C)...")
    today = datetime.date.today()

    # Precompute earliest downgrade per gene
    downgrade_info = {}
    for _, r in cg[cg["classification"].isin(DOWNGRADED_CLASSIFICATIONS)].iterrows():
        sym = r["gene"]
        c_date = r["classification_date"]
        c_class = r["classification"]
        c_gcep = r.get("gcep") or "ClinGen GCEP"
        c_url = r.get("report_url") or ""
        c_dis = r.get("disease") or "Genetic Disease"
        c_mid = r.get("mondo_id") or ""

        if sym not in downgrade_info:
            downgrade_info[sym] = {
                "date": c_date,
                "class": c_class,
                "gcep": c_gcep,
                "url": c_url,
                "disease": c_dis,
                "mondo_id": c_mid,
            }
        else:
            # Prioritize Refuted/Disputed over Limited
            if c_class in SEVERITY_CLASS_A and downgrade_info[sym]["class"] not in SEVERITY_CLASS_A:
                downgrade_info[sym] = {
                    "date": c_date or downgrade_info[sym]["date"],
                    "class": c_class,
                    "gcep": c_gcep,
                    "url": c_url,
                    "disease": c_dis,
                    "mondo_id": c_mid,
                }

    findings = []
    # Class A & B: Variants in downgraded genes
    plp_variants_df = pd.read_sql_query("SELECT * FROM clinvar_variants WHERE is_plp = 1", conn)
    if not plp_variants_df.empty:
        for _, vrow in plp_variants_df.iterrows():
            sym = vrow["gene_symbol"]
            if sym in downgrade_info:
                d_info = downgrade_info[sym]
                is_class_a = d_info["class"] in SEVERITY_CLASS_A
                audit_cls = "Class A" if is_class_a else "Class B"
                sev = "CRITICAL" if is_class_a else "WARNING"

                dt_diff = 0
                if d_info["date"]:
                    try:
                        d_obj = datetime.datetime.strptime(str(d_info["date"])[:10], "%Y-%m-%d").date()
                        dt_diff = max(0, (today - d_obj).days)
                    except Exception:
                        dt_diff = 365

                v_id = str(vrow["variation_id"])
                v_name = vrow["variant_name"]
                stars_val = vrow.get("stars", 0)
                n_sub = vrow.get("num_submitters", 1)
                t_stat = f"{vrow['clinical_significance']} ({stars_val} stars)"

                # Check population frequency
                cat_entry = GNOMAD_VARIANT_CATALOG.get(int(vrow["variation_id"]))
                af_val = cat_entry["af_global"] if cat_entry else 0.0
                mcaf_val = 0.0001
                pop_refuted = 1 if (cat_entry and af_val > mcaf_val) else 0

                evid = (
                    f"ClinGen GCEP ({d_info['gcep']}) classified gene '{sym}' as '{d_info['class']}' on {d_info['date']}. "
                    f"Downstream ClinVar database continues to report this variant as {vrow['clinical_significance']} "
                    f"with {stars_val} stars."
                )

                findings.append({
                    "audit_class": audit_cls,
                    "severity": sev,
                    "gene_symbol": sym,
                    "disease_id": d_info["mondo_id"] or d_info["disease"],
                    "disease_name": d_info["disease"],
                    "clingen_classification": d_info["class"],
                    "clingen_date": d_info["date"],
                    "clingen_gcep": d_info["gcep"],
                    "clingen_report_url": d_info["url"],
                    "target_type": "ClinVar Variant",
                    "target_id": v_id,
                    "target_name": v_name,
                    "target_status": t_stat,
                    "evidence_summary": evid,
                    "submitters_or_panel_info": f"{n_sub} submitter(s)",
                    "last_evaluated_date": vrow.get("last_evaluated"),
                    "days_unresolved": dt_diff,
                    "gnomad_af": af_val,
                    "mcaf_threshold": mcaf_val,
                    "population_refuted_flag": pop_refuted,
                })

    # Class C: PanelApp Green in Disputed/Refuted genes
    green_panels_df = pd.read_sql_query("SELECT * FROM panelapp_entries WHERE LOWER(confidence) = 'green'", conn)
    if not green_panels_df.empty:
        for _, prow in green_panels_df.iterrows():
            sym = prow["gene_symbol"]
            if sym in downgrade_info:
                d_info = downgrade_info[sym]
                is_class_c = d_info["class"] in SEVERITY_CLASS_A
                audit_cls = "Class C" if is_class_c else "Class B"
                sev = "CRITICAL" if is_class_c else "WARNING"

                dt_diff = 0
                if d_info["date"]:
                    try:
                        d_obj = datetime.datetime.strptime(str(d_info["date"])[:10], "%Y-%m-%d").date()
                        dt_diff = max(0, (today - d_obj).days)
                    except Exception:
                        dt_diff = 365

                p_id = str(prow.get("panel_id") or "PanelApp")
                p_name = prow.get("panel_name") or f"Panel {p_id}"
                p_stat = "Green (Diagnostic-Grade)"

                evid = (
                    f"ClinGen GCEP ({d_info['gcep']}) refuted/disputed gene validity as '{d_info['class']}'. "
                    f"However, Genomics England PanelApp still lists '{sym}' as diagnostic-grade Green in {p_name}."
                )

                findings.append({
                    "audit_class": audit_cls,
                    "severity": sev,
                    "gene_symbol": sym,
                    "disease_id": d_info["mondo_id"] or d_info["disease"],
                    "disease_name": d_info["disease"],
                    "clingen_classification": d_info["class"],
                    "clingen_date": d_info["date"],
                    "clingen_gcep": d_info["gcep"],
                    "clingen_report_url": d_info["url"],
                    "target_type": "PanelApp Panel",
                    "target_id": p_id,
                    "target_name": p_name,
                    "target_status": p_stat,
                    "evidence_summary": evid,
                    "submitters_or_panel_info": f"Panel version {prow.get('panel_version', '1.0')}",
                    "last_evaluated_date": d_info["date"],
                    "days_unresolved": dt_diff,
                    "gnomad_af": 0.0,
                    "mcaf_threshold": 0.0001,
                    "population_refuted_flag": 0,
                })

    if findings:
        findings_df = pd.DataFrame(findings)
        findings_df.to_sql("audit_findings", conn, if_exists="append", index=False)
        print(f"         Materialized {len(findings_df):,} audit inconsistency findings across Classes A, B, and C.")

    # 6. Legacy Views & B-Tree Indexing
    print("   [6/6] Creating compatibility views and B-tree performance indexes...")
    conn.executescript(
        """
        -- Legacy table aliases
        CREATE VIEW IF NOT EXISTS gene AS SELECT * FROM genes;
        CREATE VIEW IF NOT EXISTS clingen_curation AS SELECT * FROM clingen_curations;
        CREATE VIEW IF NOT EXISTS clinvar_variant AS SELECT * FROM clinvar_variants;
        CREATE VIEW IF NOT EXISTS panelapp_entry AS SELECT * FROM panelapp_entries;

        -- Legacy cascade views
        CREATE VIEW IF NOT EXISTS v_cascade AS
        SELECT  g.gene_symbol,
                c.classification,
                MIN(c.classification_date)                        AS assertion_date,
                COUNT(DISTINCT v.variation_id)                    AS plp_variants,
                SUM(CASE WHEN v.stars >= 2 THEN 1 ELSE 0 END)     AS plp_2star,
                SUM(CASE WHEN v.stars >= 3 THEN 1 ELSE 0 END)     AS plp_3star,
                COUNT(DISTINCT CASE WHEN p.confidence = 'green'
                                    THEN p.panel_id END)          AS green_panels
        FROM       genes            g
        JOIN       clingen_curations c ON c.gene_symbol = g.gene_symbol
        LEFT JOIN  clinvar_variants  v ON v.gene_symbol = g.gene_symbol AND v.is_plp = 1
        LEFT JOIN  panelapp_entries   p ON p.gene_symbol = g.gene_symbol
        WHERE      g.is_downgraded = 1
        GROUP BY   g.gene_symbol, c.classification;

        CREATE VIEW IF NOT EXISTS v_inconsistent AS
        SELECT * FROM v_cascade
        WHERE  green_panels > 0 AND plp_variants > 0
        ORDER BY plp_2star DESC, plp_variants DESC;

        CREATE VIEW IF NOT EXISTS v_inconsistent_genes AS
        SELECT 
            af.gene_symbol,
            g.hgnc_id,
            af.disease_id,
            af.disease_name,
            af.clingen_classification,
            af.clingen_date,
            af.clingen_gcep,
            COUNT(DISTINCT CASE WHEN af.target_type = 'ClinVar Variant' THEN af.target_id END) as total_plp_variants,
            COUNT(DISTINCT CASE WHEN af.target_type = 'PanelApp Panel' THEN af.target_id END) as total_panels,
            MAX(af.days_unresolved) as max_days_unresolved
        FROM audit_findings af
        LEFT JOIN genes g ON af.gene_symbol = g.gene_symbol
        GROUP BY af.gene_symbol, af.disease_id;

        CREATE VIEW IF NOT EXISTS v_gene_cascade_summary AS
        SELECT * FROM v_cascade;

        -- B-Tree Performance Indexes
        CREATE INDEX IF NOT EXISTS idx_full_genes_down ON genes(is_downgraded);
        CREATE INDEX IF NOT EXISTS idx_full_cur_gene ON clingen_curations(gene_symbol);
        CREATE INDEX IF NOT EXISTS idx_full_cur_class ON clingen_curations(classification);
        CREATE INDEX IF NOT EXISTS idx_full_cur_disease ON clingen_curations(disease_id);
        CREATE INDEX IF NOT EXISTS idx_full_var_gene ON clinvar_variants(gene_symbol);
        CREATE INDEX IF NOT EXISTS idx_full_var_stars ON clinvar_variants(stars);
        CREATE INDEX IF NOT EXISTS idx_full_var_plp ON clinvar_variants(is_plp);
        CREATE INDEX IF NOT EXISTS idx_full_pan_gene ON panelapp_entries(gene_symbol);
        CREATE INDEX IF NOT EXISTS idx_full_pan_conf ON panelapp_entries(confidence);
        CREATE INDEX IF NOT EXISTS idx_full_find_gene ON audit_findings(gene_symbol);
        CREATE INDEX IF NOT EXISTS idx_full_find_class ON audit_findings(audit_class);
        CREATE INDEX IF NOT EXISTS idx_full_find_gene_dis ON audit_findings(gene_symbol, disease_id);
        """
    )

    conn.commit()
    conn.close()

    # 7. Enrich Population Genetics & Multi-Ancestry
    print("   Enriching population genetics & multi-ancestry metrics...")
    enrich_database_population_genetics(db_path)

    # 8. Sync to curacascade/data/out/curacascade.db if requested
    if sync_compat and COMPAT_FULL_DB != db_path:
        try:
            import shutil
            COMPAT_FULL_DB.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(db_path, COMPAT_FULL_DB)
            print(f"   Synchronized compatibility database: {COMPAT_FULL_DB}")
        except Exception as e:
            print(f"   Note: compatibility copy skipped: {e}")

    elapsed = time.time() - start_time
    print(f"✅ [ETL Complete] Built full-scale warehouse in {elapsed:.1f}s ({db_path.stat().st_size / 1e6:.1f} MB)")

    # Summary Report
    with sqlite3.connect(str(db_path)) as c_verify:
        cur = c_verify.cursor()
        counts = {
            "total_genes": cur.execute("SELECT COUNT(*) FROM genes").fetchone()[0],
            "downgraded_genes": cur.execute("SELECT COUNT(*) FROM genes WHERE is_downgraded = 1").fetchone()[0],
            "total_curations": cur.execute("SELECT COUNT(*) FROM clingen_curations").fetchone()[0],
            "total_variants": cur.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0],
            "total_panels": cur.execute("SELECT COUNT(*) FROM panelapp_entries").fetchone()[0],
            "total_findings": cur.execute("SELECT COUNT(*) FROM audit_findings").fetchone()[0],
            "class_a_findings": cur.execute("SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class A'").fetchone()[0],
            "class_b_findings": cur.execute("SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class B'").fetchone()[0],
            "class_c_findings": cur.execute("SELECT COUNT(*) FROM audit_findings WHERE audit_class = 'Class C'").fetchone()[0],
            "elapsed_seconds": round(elapsed, 1),
            "db_path": str(db_path),
        }
        print(f"   📊 Statistics: {counts['total_genes']:,} Genes | {counts['total_findings']:,} Findings ({counts['class_a_findings']:,} Class A, {counts['class_b_findings']:,} Class B, {counts['class_c_findings']:,} Class C)")

    return counts


if __name__ == "__main__":
    build_full_warehouse()
