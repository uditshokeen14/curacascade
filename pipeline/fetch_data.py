"""CuraCascade — Real Data Fetchers and Stream Ingestion Parsers.

Implements network-resilient ingestion tools for:
1. ClinGen Gene-Disease Validity CSV:
   - Direct download from public search URL or local fallback.
   - Preamble skipping and column normalization.
2. Genomics England PanelApp API Client:
   - Queries PanelApp REST API v1 for target clinical panels.
   - Normalizes confidence tiers ('green', 'amber', 'red') and inheritance.
3. ClinVar Streamed Ingestion:
   - Streamed, memory-bounded extraction from variant_summary.txt.gz.
   - Column resolution strictly by header name (never positional index).
   - In-flight filtering matching only target curated genes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Union

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd
import requests

import config as C
from database.connection import (
    ensure_database_ready,
    execute_many,
)
from pipeline.audit import run_audit

# ==============================================================================
# 1. ClinGen Gene-Disease Validity Ingestion
# ==============================================================================

def fetch_clingen_data(
    source: Optional[Union[str, Path]] = None,
    save_raw: bool = True,
) -> pd.DataFrame:
    """Fetch and parse ClinGen Gene-Disease Validity assertions.

    Skips metadata preamble, handles delimiters, and standardizes columns:
    gene_symbol, hgnc_id, disease_name, disease_id, moi, classification,
    classification_date, gcep, report_url.

    Args:
        source: Custom URL or local file path. If None, uses default ClinGen URL.
        save_raw: If True, caches raw content to data/raw/clingen_validity_raw.csv.

    Returns:
        Cleaned Pandas DataFrame of ClinGen curations.
    """
    raw_text: str = ""

    # Check local fallback in existing curacascade/data/interim if no source specified
    fallback_interim = C.BASE_DIR / "curacascade" / "data" / "interim" / "clingen_validity.csv.gz"
    
    if source is not None:
        src_path = Path(source)
        if src_path.exists():
            print(f"Reading ClinGen data from local file: {src_path}")
            if str(src_path).endswith(".gz"):
                with gzip.open(src_path, "rt", encoding="utf-8", errors="replace") as f:
                    raw_text = f.read()
            else:
                raw_text = src_path.read_text(encoding="utf-8", errors="replace")
        else:
            # Treat as URL
            url = str(source)
            print(f"Downloading ClinGen validity data from URL: {url}")
            resp = requests.get(url, headers=C.REQUEST_HEADERS, timeout=C.REQUEST_TIMEOUT)
            resp.raise_for_status()
            raw_text = resp.text
    elif fallback_interim.exists():
        print(f"Reading ClinGen data from cached project file: {fallback_interim}")
        with gzip.open(fallback_interim, "rt", encoding="utf-8", errors="replace") as f:
            raw_text = f.read()
    else:
        url = C.CLINGEN_VALIDITY_URL
        print(f"Fetching ClinGen validity snapshot from: {url}")
        try:
            resp = requests.get(url, headers=C.REQUEST_HEADERS, timeout=C.REQUEST_TIMEOUT)
            resp.raise_for_status()
            raw_text = resp.text
        except Exception as err:
            print(f"ClinGen download failed ({err}). Checking data/raw/ ...")
            cached_raw = C.RAW_DATA_DIR / "clingen_validity_raw.csv"
            if cached_raw.exists():
                print(f"Using cached raw file: {cached_raw}")
                raw_text = cached_raw.read_text(encoding="utf-8", errors="replace")
            else:
                raise RuntimeError(
                    f"Unable to fetch ClinGen validity data from {url} and no local cache found. "
                    "Run 'python -m pipeline.seed_sample_data' for immediate offline use."
                ) from err

    if save_raw and raw_text:
        raw_out = C.RAW_DATA_DIR / "clingen_validity_raw.csv"
        raw_out.write_text(raw_text, encoding="utf-8")

    # Parse ClinGen CSV: Skip preamble lines until the header is encountered
    lines = raw_text.splitlines()
    header_idx = -1
    for i, line in enumerate(lines):
        # Look for typical ClinGen header columns
        lower_line = line.lower()
        if "gene symbol" in lower_line or "gene" in lower_line and "classification" in lower_line:
            header_idx = i
            break

    if header_idx == -1:
        # Fallback to standard line 0 or check if already clean CSV
        header_idx = 0

    table_data = "\n".join(lines[header_idx:])
    df = pd.read_csv(io.StringIO(table_data), dtype=str)

    # Standardize column naming
    col_map = {
        "GENE SYMBOL": "gene_symbol",
        "Gene Symbol": "gene_symbol",
        "gene": "gene_symbol",
        "HGNC ID": "hgnc_id",
        "hgnc_id": "hgnc_id",
        "DISEASE LABEL": "disease_name",
        "Disease Label": "disease_name",
        "disease": "disease_name",
        "DISEASE ID (MONDO)": "disease_id",
        "Disease ID": "disease_id",
        "mondo_id": "disease_id",
        "MOI": "moi",
        "moi": "moi",
        "CLASSIFICATION": "classification",
        "Classification": "classification",
        "classification": "classification",
        "CLASSIFICATION DATE": "classification_date",
        "Classification Date": "classification_date",
        "classification_date": "classification_date",
        "GCEP": "gcep",
        "gcep": "gcep",
        "CLASSIFICATION REPORT URL": "report_url",
        "Classification Report URL": "report_url",
        "report_url": "report_url",
        "SOP": "sop",
        "sop": "sop",
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Drop any decorative separator rows (e.g. +++++++++++)
    if "gene_symbol" in df.columns:
        df = df[df["gene_symbol"].notna()]
        df = df[~df["gene_symbol"].str.startswith("+")]
        df["gene_symbol"] = df["gene_symbol"].str.strip().str.upper()

    # Format classification date to YYYY-MM-DD
    if "classification_date" in df.columns:
        df["classification_date"] = pd.to_datetime(
            df["classification_date"], errors="coerce"
        ).dt.strftime("%Y-%m-%d")

    # Clean classification strings
    if "classification" in df.columns:
        df["classification"] = df["classification"].str.strip()

    # Default missing fields if not present
    for req_col in ("hgnc_id", "disease_name", "disease_id", "moi", "gcep", "report_url", "sop"):
        if req_col not in df.columns:
            df[req_col] = None

    print(f"Parsed {len(df):,} ClinGen gene-disease validity curations.")
    return df

# ==============================================================================
# 2. Genomics England PanelApp API Client
# ==============================================================================

def fetch_panelapp_panel(
    panel_id: int = 134,
    save_raw: bool = True,
) -> pd.DataFrame:
    """Fetch gene membership from Genomics England PanelApp REST API v1.

    Args:
        panel_id: Numerical PanelApp identifier (e.g. 134 for Adult HCM, 137 for LQTS, 55 for HBOC).
        save_raw: If True, writes raw JSON to data/raw/.

    Returns:
        DataFrame containing gene_symbol, panel_id, panel_name, panel_version,
        confidence, mode_of_inheritance, phenotypes.
    """
    url = f"{C.PANELAPP_BASE_URL}/panels/{panel_id}/"
    print(f"Querying PanelApp API: {url} ...")

    try:
        resp = requests.get(url, headers=C.REQUEST_HEADERS, timeout=C.REQUEST_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
    except Exception as err:
        print(f"PanelApp API query failed for panel {panel_id}: {err}")
        # Check cache or fallback
        raw_cache = C.RAW_DATA_DIR / f"panelapp_panel_{panel_id}.json"
        if raw_cache.exists():
            print(f"Loading cached PanelApp JSON from: {raw_cache}")
            data = json.loads(raw_cache.read_text(encoding="utf-8"))
        else:
            print("No cached JSON found for panel. Returning empty DataFrame.")
            return pd.DataFrame()

    if save_raw and data:
        raw_out = C.RAW_DATA_DIR / f"panelapp_panel_{panel_id}.json"
        raw_out.write_text(json.dumps(data, indent=2), encoding="utf-8")

    panel_name = data.get("name", f"Panel {panel_id}")
    panel_version = str(data.get("version", "1.0"))
    genes_data = data.get("genes", [])

    rows: List[Dict[str, Any]] = []
    for g_entry in genes_data:
        gene_info = g_entry.get("gene_data", {})
        symbol = (gene_info.get("gene_symbol") or "").strip().upper()
        if not symbol:
            continue

        conf_level = str(g_entry.get("confidence_level", "")).lower()
        # Map numeric or textual confidence:
        # 3 -> green, 2 -> amber, 1 -> red, 0 -> red
        if conf_level == "3" or "green" in conf_level:
            confidence = "green"
        elif conf_level == "2" or "amber" in conf_level:
            confidence = "amber"
        else:
            confidence = "red"

        moi = g_entry.get("mode_of_inheritance", "")
        phenotypes_list = g_entry.get("phenotypes", [])
        phenotypes = "; ".join(phenotypes_list) if isinstance(phenotypes_list, list) else str(phenotypes_list)

        rows.append({
            "gene_symbol": symbol,
            "panel_id": panel_id,
            "panel_name": panel_name,
            "panel_version": panel_version,
            "confidence": confidence,
            "mode_of_inheritance": moi,
            "phenotypes": phenotypes,
        })

    df = pd.DataFrame(rows)
    print(f"Retrieved {len(df)} genes for panel '{panel_name}' (ID: {panel_id}).")
    return df


def fetch_standard_panels(panel_ids: Optional[List[int]] = None) -> pd.DataFrame:
    """Fetch multiple standard clinical panels and combine them.

    Args:
        panel_ids: List of panel IDs. Defaults to [134, 137, 143, 55, 234, 189].

    Returns:
        Combined DataFrame of panel memberships.
    """
    if panel_ids is None:
        panel_ids = [134, 137, 143, 55, 234, 189]

    frames = []
    for pid in panel_ids:
        try:
            pdf = fetch_panelapp_panel(pid)
            if not pdf.empty:
                frames.append(pdf)
            time.sleep(0.5)  # Courteous pacing
        except Exception as e:
            print(f"Skipping panel {pid} due to error: {e}")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).drop_duplicates(subset=["gene_symbol", "panel_id"])

# ==============================================================================
# 3. ClinVar Streamed / Filtered Ingestion
# ==============================================================================

def is_pathogenic_or_likely_pathogenic(significance: str) -> bool:
    """Determine if a ClinVar clinical significance assertion is P or LP.

    Explicitly excludes 'Conflicting classifications of pathogenicity' and 'VUS'.

    Args:
        significance: ClinVar ClinicalSignificance string.

    Returns:
        True if assertion is Pathogenic or Likely Pathogenic, False otherwise.
    """
    if not significance or pd.isna(significance):
        return False
    s_lower = str(significance).lower().strip()
    if "conflicting" in s_lower:
        return False
    for kw in C.PLP_KEYWORDS:
        if kw in s_lower:
            return True
    return False


def map_review_status_to_stars(review_status: str) -> int:
    """Map a ClinVar review status string to 0-4 star rating.

    Args:
        review_status: ClinVar ReviewStatus string.

    Returns:
        Integer star rating from 0 to 4.
    """
    if not review_status or pd.isna(review_status):
        return 0
    norm = str(review_status).lower().strip()
    return C.REVIEW_STARS.get(norm, 0)


def stream_filter_clinvar(
    target_genes: Set[str],
    source_path_or_url: Optional[Union[str, Path]] = None,
    max_records: Optional[int] = None,
) -> pd.DataFrame:
    """Stream-filter ClinVar variant_summary keeping only variants in target genes.

    Resolves columns strictly by header name to guard against historical layout drift.

    Args:
        target_genes: Set of uppercase gene symbols to extract.
        source_path_or_url: Local file path or URL to variant_summary.txt.gz.
        max_records: Optional cap on retained variants.

    Returns:
        DataFrame of filtered ClinVar assertions.
    """
    target_set = {g.strip().upper() for g in target_genes if g}
    print(f"Stream-filtering ClinVar for {len(target_set)} target genes ...")

    # Local interim fallback in existing curacascade/data/interim
    fallback_clinvar = C.BASE_DIR / "curacascade" / "data" / "interim" / "clinvar_latest.csv.gz"

    if source_path_or_url is None and fallback_clinvar.exists():
        print(f"Using pre-filtered project ClinVar archive: {fallback_clinvar}")
        cv_df = pd.read_csv(fallback_clinvar)
        # Rename columns to standard names
        rename_dict = {
            "VariationID": "variation_id",
            "GeneSymbol": "gene_symbol",
            "Name": "variant_name",
            "Type": "variant_type",
            "ClinicalSignificance": "clinical_significance",
            "ReviewStatus": "review_status",
            "NumberSubmitters": "num_submitters",
            "LastEvaluated": "last_evaluated",
            "PhenotypeList": "phenotype_list",
            "gene": "gene_symbol",
        }
        cv_df = cv_df.rename(columns={k: v for k, v in rename_dict.items() if k in cv_df.columns})
        if "gene_symbol" in cv_df.columns:
            cv_df = cv_df[cv_df["gene_symbol"].str.upper().isin(target_set)].copy()
            cv_df["is_plp"] = cv_df["clinical_significance"].apply(is_pathogenic_or_likely_pathogenic).astype(int)
            if "stars" not in cv_df.columns and "review_status" in cv_df.columns:
                cv_df["stars"] = cv_df["review_status"].apply(map_review_status_to_stars)
            return cv_df

    # Stream from local or remote variant_summary.txt.gz
    path = Path(source_path_or_url) if source_path_or_url else None
    if path and path.exists():
        file_obj = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        url = str(source_path_or_url) if source_path_or_url else C.CLINVAR_SUMMARY_URL
        print(f"Streaming ClinVar over HTTP from: {url}")
        resp = requests.get(url, stream=True, headers=C.REQUEST_HEADERS, timeout=60)
        resp.raise_for_status()
        file_obj = gzip.GzipFile(fileobj=resp.raw, mode="rb")
        file_obj = io.TextIOWrapper(file_obj, encoding="utf-8", errors="replace")

    extracted_variants: List[Dict[str, Any]] = []

    try:
        reader = csv.DictReader(file_obj, delimiter="\t")
        for row in reader:
            gene = (row.get("GeneSymbol") or row.get("gene") or "").strip().upper()
            if gene not in target_set:
                continue

            var_id_str = row.get("VariationID") or row.get("variation_id") or "0"
            try:
                var_id = int(var_id_str)
            except ValueError:
                continue

            sig = row.get("ClinicalSignificance", "")
            rev_stat = row.get("ReviewStatus", "")
            last_eval = row.get("LastEvaluated", "")
            if last_eval == "-" or not last_eval:
                last_eval = None

            num_sub_str = row.get("NumberSubmitters", "1")
            try:
                num_sub = int(num_sub_str) if num_sub_str and num_sub_str != "-" else 1
            except ValueError:
                num_sub = 1

            extracted_variants.append({
                "variation_id": var_id,
                "gene_symbol": gene,
                "variant_name": row.get("Name", f"Variation {var_id}"),
                "variant_type": row.get("Type", "unknown"),
                "clinical_significance": sig,
                "is_plp": 1 if is_pathogenic_or_likely_pathogenic(sig) else 0,
                "review_status": rev_stat,
                "stars": map_review_status_to_stars(rev_stat),
                "num_submitters": num_sub,
                "submitter_names": row.get("SubmitterNames", None),
                "last_evaluated": last_eval,
                "phenotype_list": row.get("PhenotypeList", ""),
            })

            if max_records and len(extracted_variants) >= max_records:
                print(f"Reached cap of {max_records} ClinVar variants.")
                break
    finally:
        file_obj.close()

    result_df = pd.DataFrame(extracted_variants)
    if not result_df.empty:
        result_df = result_df.drop_duplicates(subset=["variation_id"])
    print(f"Extracted {len(result_df):,} ClinVar assertions for {len(target_set)} target genes.")
    return result_df

# ==============================================================================
# 4. Relational Ingestion Manager
# ==============================================================================

def ingest_all_into_db(
    clingen_df: pd.DataFrame,
    panelapp_df: Optional[pd.DataFrame] = None,
    clinvar_df: Optional[pd.DataFrame] = None,
    db_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Populate relational tables from parsed DataFrames and run the audit.

    Args:
        clingen_df: ClinGen curations DataFrame.
        panelapp_df: Optional PanelApp memberships DataFrame.
        clinvar_df: Optional ClinVar variants DataFrame.
        db_path: Target database path.

    Returns:
        Audit report metrics dictionary.
    """
    ensure_database_ready(db_path=db_path, drop_existing=True)

    # 1. Insert Unique Diseases
    print("Normalizing diseases ...")
    diseases_dict: Dict[str, str] = {}
    if "disease_id" in clingen_df.columns and "disease_name" in clingen_df.columns:
        for _, r in clingen_df.iterrows():
            did = str(r["disease_id"]).strip() if pd.notna(r["disease_id"]) else None
            dname = str(r["disease_name"]).strip() if pd.notna(r["disease_name"]) else "Unknown Disease"
            if did and did not in diseases_dict:
                diseases_dict[did] = dname

    disease_rows = [(did, dname, did if "MONDO" in did else None, None, None) for did, dname in diseases_dict.items()]
    if disease_rows:
        execute_many(
            "INSERT OR IGNORE INTO diseases (disease_id, disease_name, mondo_id, omim_id, inheritance_mode) VALUES (?, ?, ?, ?, ?);",
            disease_rows,
            db_path=db_path,
        )

    # 2. Insert Genes
    print("Populating genes table ...")
    downgraded_genes = set(clingen_df.loc[clingen_df["classification"].isin(C.DOWNGRADED_CLASSIFICATIONS), "gene_symbol"])
    
    unique_genes = clingen_df.groupby("gene_symbol", as_index=False).agg({
        "hgnc_id": "first",
    })

    gene_rows = []
    for _, gr in unique_genes.iterrows():
        sym = gr["gene_symbol"]
        is_down = 1 if sym in downgraded_genes else 0
        # Determine highest severity for this gene
        gene_curations = clingen_df[clingen_df["gene_symbol"] == sym]["classification"].tolist()
        highest_sev = "None"
        for sev in ("Refuted", "Disputed", "Limited", "Moderate", "Strong", "Definitive"):
            if sev in gene_curations:
                highest_sev = sev
                break
        gene_rows.append((sym, gr.get("hgnc_id"), None, None, None, is_down, highest_sev))

    execute_many(
        "INSERT OR IGNORE INTO genes (gene_symbol, hgnc_id, entrez_id, gene_name, chromosome, is_downgraded, highest_downgrade_severity) VALUES (?, ?, ?, ?, ?, ?, ?);",
        gene_rows,
        db_path=db_path,
    )

    # 3. Insert ClinGen Curations
    print("Populating clingen_curations ...")
    curation_rows = []
    for _, cr in clingen_df.iterrows():
        curation_rows.append((
            cr["gene_symbol"],
            cr.get("disease_id") or "UNKNOWN",
            cr.get("moi"),
            cr["classification"],
            cr.get("classification_date"),
            cr.get("gcep"),
            cr.get("sop"),
            cr.get("report_url"),
            None,
        ))
    execute_many(
        "INSERT INTO clingen_curations (gene_symbol, disease_id, moi, classification, classification_date, gcep, sop, report_url, pmids) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
        curation_rows,
        db_path=db_path,
    )

    # 4. Insert ClinVar Variants
    if clinvar_df is not None and not clinvar_df.empty:
        print(f"Populating clinvar_variants ({len(clinvar_df):,} records) ...")
        valid_genes = set(unique_genes["gene_symbol"])
        var_rows = []
        for _, vr in clinvar_df.iterrows():
            sym = vr["gene_symbol"]
            if sym not in valid_genes:
                continue
            var_rows.append((
                int(vr["variation_id"]),
                sym,
                vr.get("disease_id"),
                vr.get("variant_name", f"Var {vr['variation_id']}"),
                vr.get("variant_type", "single nucleotide variant"),
                vr.get("clinical_significance", "Unknown"),
                int(vr.get("is_plp", 0)),
                vr.get("review_status", "no assertion criteria provided"),
                int(vr.get("stars", 0)),
                int(vr.get("num_submitters", 1)),
                vr.get("submitter_names"),
                vr.get("last_evaluated"),
                None,
                vr.get("phenotype_list"),
            ))
        execute_many(
            """
            INSERT OR IGNORE INTO clinvar_variants (
                variation_id, gene_symbol, disease_id, variant_name, variant_type, clinical_significance,
                is_plp, review_status, stars, num_submitters, submitter_names, last_evaluated, rcv_accessions, phenotype_list
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """,
            var_rows,
            db_path=db_path,
        )

    # 5. Insert PanelApp Entries
    if panelapp_df is not None and not panelapp_df.empty:
        print(f"Populating panelapp_entries ({len(panelapp_df):,} records) ...")
        valid_genes = set(unique_genes["gene_symbol"])
        pan_rows = []
        for _, pr in panelapp_df.iterrows():
            sym = pr["gene_symbol"]
            if sym not in valid_genes:
                continue
            pan_rows.append((
                sym,
                None,
                int(pr["panel_id"]),
                pr.get("panel_name", f"Panel {pr['panel_id']}"),
                pr.get("panel_version", "1.0"),
                str(pr.get("confidence", "amber")).lower(),
                pr.get("mode_of_inheritance"),
                pr.get("phenotypes"),
            ))
        execute_many(
            """
            INSERT INTO panelapp_entries (
                gene_symbol, disease_id, panel_id, panel_name, panel_version, confidence, mode_of_inheritance, phenotypes
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?);
            """,
            pan_rows,
            db_path=db_path,
        )

    print("Running audit on newly ingested data ...")
    audit_results = run_audit(db_path=db_path)
    return audit_results


def main():
    """CLI entry point for data fetching and ingestion."""
    parser = argparse.ArgumentParser(description="Fetch and ingest real genomic databases into CuraCascade.")
    parser.add_argument("--clingen", action="store_true", help="Fetch ClinGen validity dataset.")
    parser.add_argument("--panelapp", action="store_true", help="Fetch PanelApp panels.")
    parser.add_argument("--clinvar", action="store_true", help="Stream-filter ClinVar variants.")
    parser.add_argument("--all", action="store_true", help="Execute complete ingestion workflow.")
    parser.add_argument("--db", type=str, default=str(C.DB_PATH), help="Target SQLite database.")
    args = parser.parse_args()

    run_all = args.all or not (args.clingen or args.panelapp or args.clinvar)

    print("=== CuraCascade Data Ingestion Pipeline ===")
    cg_df = fetch_clingen_data()

    pa_df = pd.DataFrame()
    if run_all or args.panelapp:
        pa_df = fetch_standard_panels()

    cv_df = pd.DataFrame()
    if run_all or args.clinvar:
        target_genes = set(cg_df["gene_symbol"])
        cv_df = stream_filter_clinvar(target_genes)

    ingest_all_into_db(cg_df, pa_df, cv_df, db_path=args.db)
    print("Ingestion complete!")


if __name__ == "__main__":
    main()
