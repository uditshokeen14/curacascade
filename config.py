"""CuraCascade — Configuration Module.

Defines file paths, database connection parameters, external API endpoints,
ClinGen curation taxonomies, ClinVar review status star mappings,
and inconsistency audit classification thresholds.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Final, List, Set

# ==============================================================================
# Path Configurations
# ==============================================================================
BASE_DIR: Final[Path] = Path(os.environ.get("CURACASCADE_HOME") or Path(__file__).resolve().parent)
DATA_DIR: Final[Path] = BASE_DIR / "data"
RAW_DATA_DIR: Final[Path] = DATA_DIR / "raw"
CACHE_DIR: Final[Path] = DATA_DIR / "cache"
PROCESSED_DATA_DIR: Final[Path] = DATA_DIR / "processed"

# Primary database file (SQLite3)
DB_PATH: Final[Path] = DATA_DIR / "curacascade.db"

# Create directories if they do not exist
for _dir in (DATA_DIR, RAW_DATA_DIR, CACHE_DIR, PROCESSED_DATA_DIR):
    _dir.mkdir(parents=True, exist_ok=True)

# ==============================================================================
# External Resource URLs & Endpoints
# ==============================================================================
# ClinGen public Gene-Disease Validity snapshot download URL
CLINGEN_VALIDITY_URL: Final[str] = (
    "https://search.clinicalgenome.org/kb/gene-validity/download"
)

# Genomics England PanelApp REST API v1
PANELAPP_BASE_URL: Final[str] = "https://panelapp.genomicsengland.co.uk/api/v1"

# NCBI ClinVar variant summary download URL
CLINVAR_SUMMARY_URL: Final[str] = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

# User-Agent header for courteous academic queries
USER_AGENT: Final[str] = (
    "CuraCascade/1.0 (Biological Databases Academic Project; Clinical Gene Downgrade Auditor)"
)
REQUEST_HEADERS: Final[Dict[str, str]] = {
    "User-Agent": USER_AGENT,
    "Accept": "application/json, text/csv, */*",
}
REQUEST_TIMEOUT: Final[int] = 30  # seconds

# ==============================================================================
# ClinGen Validity Taxonomy
# ==============================================================================
DOWNGRADED_CLASSIFICATIONS: Final[Set[str]] = {"Limited", "Disputed", "Refuted"}
SUPPORTED_CLASSIFICATIONS: Final[Set[str]] = {"Definitive", "Strong", "Moderate"}

SEVERITY_ORDER: Final[Dict[str, int]] = {
    "Refuted": 1,
    "Disputed": 2,
    "Limited": 3,
    "Moderate": 4,
    "Strong": 5,
    "Definitive": 6,
    "No Known Disease Relationship": 7,
}

# ==============================================================================
# ClinVar Review Status -> Star Rating
# ==============================================================================
REVIEW_STARS: Final[Dict[str, int]] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
    "no classifications from unflagged records": 0,
}

# Keywords indicating Pathogenic or Likely Pathogenic assertions in ClinVar
PLP_KEYWORDS: Final[List[str]] = [
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "pathogenic; risk factor",
    "likely pathogenic; risk factor",
]

# ==============================================================================
# Audit Inconsistency Classes
# ==============================================================================
AUDIT_CLASS_A: Final[str] = "Class A"  # P/LP in Disputed or Refuted genes
AUDIT_CLASS_B: Final[str] = "Class B"  # P/LP in Limited-evidence genes
AUDIT_CLASS_C: Final[str] = "Class C"  # PanelApp Green in Disputed or Refuted genes

CLASS_METADATA: Final[Dict[str, Dict[str, str]]] = {
    AUDIT_CLASS_A: {
        "title": "ClinVar P/LP in Disputed/Refuted Gene",
        "severity": "CRITICAL",
        "description": (
            "ClinVar reports Pathogenic or Likely Pathogenic variants for a gene whose "
            "clinical validity has been formally Disputed or Refuted by ClinGen GCEPs."
        ),
        "clinical_risk": "High risk of false-positive molecular diagnosis in a non-causal gene.",
        "badge_color": "#dc2626",  # Red
    },
    AUDIT_CLASS_B: {
        "title": "ClinVar P/LP in Limited-Evidence Gene",
        "severity": "HIGH",
        "description": (
            "ClinVar reports Pathogenic or Likely Pathogenic variants in a gene with only "
            "Limited evidence in ClinGen. According to ACMG standards, variants in Limited genes "
            "should not exceed VUS classification."
        ),
        "clinical_risk": "Unwarranted diagnostic weight given to variants in genes lacking confirmed causality.",
        "badge_color": "#ea580c",  # Orange
    },
    AUDIT_CLASS_C: {
        "title": "PanelApp Green in Disputed/Refuted Gene",
        "severity": "CRITICAL",
        "description": (
            "Genomics England PanelApp lists the gene as 'Green' (diagnostic-grade / routine reporting) "
            "for a disease where ClinGen has Disputed or Refuted gene-disease causality."
        ),
        "clinical_risk": "Diagnostic panels routinely screening and reporting genes refuted by expert consortia.",
        "badge_color": "#7c3aed",  # Purple
    },
}
