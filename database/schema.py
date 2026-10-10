"""CuraCascade — Relational Database Schema Definition.

Implements clean 3NF relational tables with strict primary keys, foreign keys,
and multi-column indexes on (gene_symbol, disease_id) for rapid join performance.
"""

from __future__ import annotations

import sqlite3
from typing import Final

DDL_SCHEMA: Final[str] = """
PRAGMA foreign_keys = ON;

-- -----------------------------------------------------------------------------
-- 1. GENES Table
-- Central entity representing curated human genes.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS genes (
    gene_symbol                TEXT PRIMARY KEY,
    hgnc_id                    TEXT,
    entrez_id                  TEXT,
    gene_name                  TEXT,
    chromosome                 TEXT,
    is_downgraded              INTEGER NOT NULL DEFAULT 0,  -- 1 if Limited, Disputed, or Refuted
    highest_downgrade_severity TEXT DEFAULT 'None',         -- 'Refuted', 'Disputed', 'Limited', or 'None'
    updated_at                 TEXT DEFAULT (datetime('now'))
);

-- -----------------------------------------------------------------------------
-- 2. DISEASES Table
-- Normalized ontology entities (e.g. MONDO / OMIM disease phenotypes).
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS diseases (
    disease_id                 TEXT PRIMARY KEY,            -- e.g. 'MONDO:0005045'
    disease_name               TEXT NOT NULL,
    mondo_id                   TEXT,
    omim_id                    TEXT,
    inheritance_mode           TEXT,
    prevalence_per_million     REAL DEFAULT 1000.0,         -- Orphanet point prevalence
    orphanet_id                TEXT,                        -- e.g. 'ORPHA:217'
    updated_at                 TEXT DEFAULT (datetime('now'))
);

-- -----------------------------------------------------------------------------
-- 3. CLINGEN_CURATIONS Table
-- Expert gene-disease clinical validity assertions from ClinGen GCEPs.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS clingen_curations (
    curation_id                INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol                TEXT NOT NULL,
    disease_id                 TEXT NOT NULL,
    moi                        TEXT,                        -- Mode of Inheritance (AD, AR, XL)
    classification             TEXT NOT NULL,               -- Definitive, Strong, Moderate, Limited, Disputed, Refuted
    classification_date        TEXT,                        -- ISO date (YYYY-MM-DD); starts latency clock
    gcep                       TEXT,                        -- Gene Curation Expert Panel
    sop                        TEXT,                        -- Standard Operating Procedure version
    report_url                 TEXT,                        -- ClinGen Evidence Report URL
    pmids                      TEXT,                        -- Associated PubMed citations
    created_at                 TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (gene_symbol)  REFERENCES genes(gene_symbol) ON DELETE CASCADE,
    FOREIGN KEY (disease_id)   REFERENCES diseases(disease_id) ON DELETE CASCADE
);

-- Indexes for ClinGen queries
CREATE INDEX IF NOT EXISTS idx_clingen_gene_disease 
    ON clingen_curations(gene_symbol, disease_id);
CREATE INDEX IF NOT EXISTS idx_clingen_classification 
    ON clingen_curations(classification);
CREATE INDEX IF NOT EXISTS idx_clingen_date 
    ON clingen_curations(classification_date);

-- -----------------------------------------------------------------------------
-- 4. CLINVAR_VARIANTS Table
-- Downstream variant assertions in ClinVar.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS clinvar_variants (
    variation_id               INTEGER PRIMARY KEY,
    gene_symbol                TEXT NOT NULL,
    disease_id                 TEXT,                        -- Linked disease or phenotype
    variant_name               TEXT NOT NULL,
    variant_type               TEXT,                        -- SNV, deletion, duplication, etc.
    clinical_significance      TEXT NOT NULL,               -- Pathogenic, Likely pathogenic, VUS, etc.
    is_plp                     INTEGER NOT NULL DEFAULT 0,  -- 1 if Pathogenic or Likely Pathogenic
    review_status              TEXT NOT NULL,
    stars                      INTEGER NOT NULL DEFAULT 0,  -- 0 to 4 stars ClinVar confidence
    num_submitters             INTEGER DEFAULT 1,
    submitter_names            TEXT,
    last_evaluated             TEXT,                        -- YYYY-MM-DD
    rcv_accessions             TEXT,
    phenotype_list             TEXT,
    gnomad_af_global           REAL DEFAULT 0.0,            -- gnomAD v4 global allele frequency
    gnomad_af_popmax           REAL DEFAULT 0.0,            -- gnomAD v4 ancestry popmax allele frequency
    gnomad_popmax_ancestry     TEXT,                        -- e.g. EAS, SAS, AFR, NFE
    mcaf_threshold             REAL DEFAULT 0.0001,         -- Maximum Credible Allele Frequency threshold
    population_refuted_flag    INTEGER DEFAULT 0,           -- 1 if gnomad_af > mcaf_threshold
    created_at                 TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (gene_symbol)  REFERENCES genes(gene_symbol) ON DELETE CASCADE,
    FOREIGN KEY (disease_id)   REFERENCES diseases(disease_id) ON DELETE SET NULL
);

-- Indexes for ClinVar queries
CREATE INDEX IF NOT EXISTS idx_clinvar_gene_disease 
    ON clinvar_variants(gene_symbol, disease_id);
CREATE INDEX IF NOT EXISTS idx_clinvar_gene 
    ON clinvar_variants(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_clinvar_is_plp 
    ON clinvar_variants(is_plp);
CREATE INDEX IF NOT EXISTS idx_clinvar_stars 
    ON clinvar_variants(stars);

-- -----------------------------------------------------------------------------
-- 5. PANELAPP_ENTRIES Table
-- Downstream diagnostic gene panel memberships from Genomics England PanelApp.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS panelapp_entries (
    entry_id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol                TEXT NOT NULL,
    disease_id                 TEXT,                        -- Associated panel disease/indication
    panel_id                   INTEGER NOT NULL,
    panel_name                 TEXT NOT NULL,
    panel_version              TEXT,
    confidence                 TEXT NOT NULL,               -- 'green' (diagnostic-grade), 'amber', 'red'
    mode_of_inheritance        TEXT,
    phenotypes                 TEXT,
    created_at                 TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (gene_symbol)  REFERENCES genes(gene_symbol) ON DELETE CASCADE,
    FOREIGN KEY (disease_id)   REFERENCES diseases(disease_id) ON DELETE SET NULL
);

-- Indexes for PanelApp queries
CREATE INDEX IF NOT EXISTS idx_panelapp_gene_disease 
    ON panelapp_entries(gene_symbol, disease_id);
CREATE INDEX IF NOT EXISTS idx_panelapp_gene 
    ON panelapp_entries(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_panelapp_confidence 
    ON panelapp_entries(confidence);

-- -----------------------------------------------------------------------------
-- 6. AUDIT_FINDINGS Table
-- Materialized findings generated by the cross-database cascade auditor.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS audit_findings (
    finding_id                 INTEGER PRIMARY KEY AUTOINCREMENT,
    audit_class                TEXT NOT NULL,               -- 'Class A', 'Class B', 'Class C'
    severity                   TEXT NOT NULL,               -- 'CRITICAL', 'HIGH', 'MEDIUM'
    gene_symbol                TEXT NOT NULL,
    disease_id                 TEXT,
    disease_name               TEXT,
    clingen_classification     TEXT NOT NULL,               -- 'Disputed', 'Refuted', 'Limited'
    clingen_date               TEXT,
    clingen_gcep               TEXT,
    clingen_report_url         TEXT,
    target_type                TEXT NOT NULL,               -- 'ClinVar Variant' or 'PanelApp Panel'
    target_id                  TEXT NOT NULL,               -- VariationID or PanelID
    target_name                TEXT NOT NULL,               -- Variant name or Panel name
    target_status              TEXT NOT NULL,               -- e.g. 'Pathogenic (2 stars)' or 'Green Diagnostic-Grade'
    evidence_summary           TEXT,
    submitters_or_panel_info   TEXT,
    last_evaluated_date        TEXT,
    days_unresolved            INTEGER,
    gnomad_af                  REAL DEFAULT 0.0,
    mcaf_threshold             REAL DEFAULT 0.0,
    population_refuted_flag    INTEGER DEFAULT 0,
    created_at                 TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (gene_symbol)  REFERENCES genes(gene_symbol) ON DELETE CASCADE,
    FOREIGN KEY (disease_id)   REFERENCES diseases(disease_id) ON DELETE SET NULL
);

-- Indexes for Audit Findings
CREATE INDEX IF NOT EXISTS idx_audit_gene_disease 
    ON audit_findings(gene_symbol, disease_id);
CREATE INDEX IF NOT EXISTS idx_audit_class 
    ON audit_findings(audit_class);
CREATE INDEX IF NOT EXISTS idx_audit_severity 
    ON audit_findings(severity);
CREATE INDEX IF NOT EXISTS idx_audit_gene 
    ON audit_findings(gene_symbol);

-- -----------------------------------------------------------------------------
-- 7. DOWNGRADE_RISK_PREDICTIONS Table (Objective 5 / Module M5)
-- Machine-learning risk forecasting for currently supported gene-disease pairs.
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS downgrade_risk_predictions (
    prediction_id              INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol                TEXT NOT NULL,
    disease_name               TEXT NOT NULL,
    clingen_status             TEXT NOT NULL,               -- 'Definitive', 'Strong', 'Moderate'
    downgrade_risk_score       REAL NOT NULL,               -- 0.00 to 1.00 calibrated probability
    risk_tier                  TEXT NOT NULL,               -- 'High Risk', 'Moderate Watch', 'Stable'
    top_risk_driver            TEXT NOT NULL,               -- Primary feature driver
    last_scored_at             TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (gene_symbol)  REFERENCES genes(gene_symbol) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_risk_predictions_gene 
    ON downgrade_risk_predictions(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_risk_predictions_tier 
    ON downgrade_risk_predictions(risk_tier);
CREATE INDEX IF NOT EXISTS idx_risk_predictions_score 
    ON downgrade_risk_predictions(downgrade_risk_score);

-- -----------------------------------------------------------------------------
-- 8. VIEWS: Analytical Cross-Layer Aggregations
-- -----------------------------------------------------------------------------
DROP VIEW IF EXISTS v_gene_cascade_summary;
CREATE VIEW v_gene_cascade_summary AS
SELECT 
    g.gene_symbol,
    g.hgnc_id,
    g.is_downgraded,
    g.highest_downgrade_severity,
    cc.classification AS clingen_classification,
    cc.classification_date AS clingen_date,
    d.disease_id,
    d.disease_name,
    COUNT(DISTINCT cv.variation_id) AS total_clinvar_variants,
    SUM(CASE WHEN cv.is_plp = 1 THEN 1 ELSE 0 END) AS plp_variants_count,
    SUM(CASE WHEN cv.is_plp = 1 AND cv.stars >= 2 THEN 1 ELSE 0 END) AS high_confidence_plp_count,
    COUNT(DISTINCT CASE WHEN pa.confidence = 'green' THEN pa.panel_id END) AS green_panels_count,
    COUNT(DISTINCT pa.panel_id) AS total_panels_count
FROM genes g
JOIN clingen_curations cc ON g.gene_symbol = cc.gene_symbol
JOIN diseases d ON cc.disease_id = d.disease_id
LEFT JOIN clinvar_variants cv ON g.gene_symbol = cv.gene_symbol
LEFT JOIN panelapp_entries pa ON g.gene_symbol = pa.gene_symbol
GROUP BY g.gene_symbol, cc.curation_id, d.disease_id;

DROP VIEW IF EXISTS v_inconsistent_genes;
CREATE VIEW v_inconsistent_genes AS
SELECT * FROM v_gene_cascade_summary
WHERE is_downgraded = 1 
  AND (plp_variants_count > 0 OR green_panels_count > 0)
ORDER BY high_confidence_plp_count DESC, plp_variants_count DESC, green_panels_count DESC;
"""


def init_db(connection: sqlite3.Connection, drop_existing: bool = False) -> None:
    """Initialize the SQLite database schema.

    Args:
        connection: Active sqlite3 connection.
        drop_existing: If True, drops all existing tables before recreating.
    """
    if drop_existing:
        drop_script = """
        PRAGMA foreign_keys = OFF;
        DROP VIEW IF EXISTS v_inconsistent_genes;
        DROP VIEW IF EXISTS v_gene_cascade_summary;
        DROP TABLE IF EXISTS downgrade_risk_predictions;
        DROP TABLE IF EXISTS audit_findings;
        DROP TABLE IF EXISTS panelapp_entries;
        DROP TABLE IF EXISTS clinvar_variants;
        DROP TABLE IF EXISTS clingen_curations;
        DROP TABLE IF EXISTS diseases;
        DROP TABLE IF EXISTS genes;
        PRAGMA foreign_keys = ON;
        """
        connection.executescript(drop_script)

    connection.executescript(DDL_SCHEMA)
    connection.commit()
