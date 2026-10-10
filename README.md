# CuraCascade: Downstream Clinical Evidence Downgrade Auditor

[![Build Status](https://github.com/udit/curacascade/actions/workflows/audit_ci.yml/badge.svg)](https://github.com/udit/curacascade/actions)
[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Database](https://img.shields.io/badge/Database-SQLite3%20%7C%20DuckDB-green.svg)](https://sqlite.org/)
[![Docker](https://img.shields.io/badge/Docker-Enabled-2496ED.svg?logo=docker)](https://www.docker.com/)
[![Dashboard](https://img.shields.io/badge/UI-Streamlit-red.svg)](https://streamlit.io/)

## Executive Abstract

**CuraCascade** is a production-grade bioinformatics database auditing framework designed to detect downstream clinical discrepancies when upstream gene–disease validity relationships are downgraded by expert consortia like **ClinGen** (*Clinical Genome Resource*). It provides robust CI/CD integration, containerized deployment, and standardized clinical dossier exports to streamline laboratory quality control workflows.

When ClinGen Gene Curation Expert Panels (GCEPs) re-appraise a gene from *Definitive* to **Disputed**, **Refuted**, or **Limited**, downstream clinical databases—specifically **NCBI ClinVar** (variant assertions) and **Genomics England PanelApp** (diagnostic panels)—frequently experience synchronization lag. CuraCascade automates cross-database relational audits to identify these critical evidence cascades.

---

## 🎯 Audit Taxonomy & Classification

| Inconsistency Class | Severity | Definition & Clinical Criterion | Clinical Risk |
| :--- | :--- | :--- | :--- |
| **Class A** | `CRITICAL` | **Pathogenic (P) or Likely Pathogenic (LP)** ClinVar assertions in **Disputed or Refuted** genes. | High probability of false-positive molecular diagnosis in a gene proven not to cause the disease. |
| **Class B** | `HIGH` | **Pathogenic (P) or Likely Pathogenic (LP)** ClinVar assertions in **Limited-evidence** genes. | Violates ACMG recommendations; variants in Limited genes should not be classified above VUS. |
| **Class C** | `CRITICAL` | **Green (diagnostic-grade / routine reporting)** PanelApp entries in **Disputed or Refuted** genes. | Routine NHS GMS diagnostic exome panels continuing to screen and report discredited genes. |

---

## 🏗️ Relational Database Architecture (3NF)

The system models clinical evidence using a 3NF relational database schema in SQLite (with DuckDB support) enforcing foreign keys and composite indexes on `(gene_symbol, disease_id)`:

```
           +---------------------+
           |        GENES        |
           |  (gene_symbol [PK]) |
           +----------+----------+
                      | 1
                      |
        +-------------+-------------+
        | N           | N           | N
+-------v-------+ +---v---------+ +-v---------------+
|    CLINGEN    | |   CLINVAR   | |    PANELAPP     |
|   CURATIONS   | |  VARIANTS   | |     ENTRIES     |
+-------+-------+ +---+---------+ +-+---------------+
        | N           | N           | N
        +-------------+-------------+
                      |
                      | 1
           +----------v----------+
           |      DISEASES       |
           |   (disease_id [PK]) |
           +---------------------+
```

### Key Relational Design Principles:
1. **Primary Keys & Foreign Keys:** Strict relational integrity between `genes`, `diseases`, `clingen_curations`, `clinvar_variants`, and `panelapp_entries`.
2. **Composite Indexes:** `(gene_symbol, disease_id)` indexed across all child tables for rapid join performance.
3. **Audit Materialization:** The `audit_findings` table stores computed discrepancies along with provenance metadata (`clingen_gcep`, `clingen_date`, `last_evaluated_date`, `days_unresolved`, `stars`, `submitter_names`).
4. **Analytical SQL Views:**
   - `v_gene_cascade_summary`: Summarizes P/LP counts and Green panel counts per downgraded gene.
   - `v_inconsistent_genes`: Materializes genes failing both downstream layers.

---

## 🚀 Quick Setup & Execution

### 1. Containerized Deployment (Docker)
The recommended deployment path is via Docker, which provides an isolated and reproducible environment.
```bash
docker build -t curacascade:latest .
docker run -p 8501:8501 curacascade:latest
```
*The dashboard will be available at http://localhost:8501.*

### 1b. Local Installation (Alternative)
```bash
pip install -r requirements.txt
```

### 2. Seed Realistic Benchmark Dataset (Immediate Local Execution)
Populates 35+ realistic, biologically curated genes (including cardiovascular re-appraisals like *MYH7*, *TNNT2*, *MYLK2*, *MYOZ2*, *KCNE2*, *SCN1B*; neurodevelopmental genes like *SHANK3*, *CHD8*, *NLGN4X*; and cancer genes like *BRCA1*, *BRCA2*, *NBN*, *RAD50*, *XRCC2*):
```bash
python -m pipeline.seed_sample_data
```

### 3. Launch Interactive Web Dashboard
```bash
streamlit run app.py
```

---

## 💻 CLI & Pipeline Commands

### Unified Headless CLI (`run.py`)
Run evidence queries, security/integrity verification, or launch the web dashboard:
```bash
# Audit specific genes for Class A critical inconsistencies:
python run.py audit --genes ACTN2,BARD1,SCN1B --class A --format table

# Export findings to JSON or CSV:
python run.py audit --genes BARD1 --format json --output findings.json
python run.py audit --class all --format csv --output all_findings.csv

# Run full database integrity & cryptographic verification suite:
python run.py verify

# Launch interactive Streamlit dashboard on port 8501:
python run.py serve
```
> **CLI Status Codes**:
> - `0`: Success / clean audit (no critical unpropagated Class A findings flagged)
> - `1`: Critical unpropagated Class A findings flagged in audit
> - `2`: Database integrity or security verification failure

### Run Audit Engine
Executes SQL joins against the current database and prints summary metrics to stdout:
```bash
python -m pipeline.audit
```

### Fetch & Stream Real Genomic Data
Fetch live data from ClinGen, query the PanelApp REST API, or stream filter NCBI ClinVar:
```bash
# Fetch real ClinGen validity snapshot:
python -m pipeline.fetch_data --clingen

# Query Genomics England PanelApp API:
python -m pipeline.fetch_data --panelapp

# Stream filter ClinVar variant_summary.txt.gz:
python -m pipeline.fetch_data --clinvar

# Full end-to-end ingestion:
python -m pipeline.fetch_data --all
```

---

## 📁 Repository Structure

```
curacascade/
├── run.py                    # Unified headless CLI interface (audit, verify, serve)
├── config.py                 # Central URLs, paths, thresholds, star mappings
├── requirements.txt          # Python runtime dependencies
├── app.py                    # Streamlit web dashboard
├── database/
│   ├── __init__.py           # Database connection helpers
│   ├── connection.py         # Thread-safe SQLite/DuckDB connection manager
│   └── schema.py             # 3NF DDL schema, tables, foreign keys, and indexes
├── pipeline/
│   ├── __init__.py           # Pipeline entrypoint exports
│   ├── seed_sample_data.py   # 35+ curated realistic genes benchmark generator
│   ├── fetch_data.py         # Real fetchers/parsers (ClinGen, PanelApp, ClinVar)
│   └── audit.py              # Cross-layer SQL joins, latency, and categorization
└── data/
    ├── raw/                  # Downloaded raw files and API cache
    ├── processed/            # Exported audit CSV reports
    └── curacascade.db        # Primary SQLite relational database
```

---

## 📊 Dashboard Features
- **KPI Summary Cards:** Total Curated Genes, Downgraded Genes, Class A, Class B, and Class C Inconsistency counts.
- **Inconsistency Audit Table:** Multi-class filtering (Class A, B, C), live gene search, severity badges, and one-click **CSV export** for clinical review.
- **Gene Detail & Timeline View:** Synchronized multi-panel inspection showing ClinGen appraisal dates, ClinVar variants (with star ratings and submitters), and PanelApp diagnostic panels.
- **Analytics View:** Distribution charts for classes, ClinVar review stars, and latency (days unresolved).

---

## 📦 Data & Software Availability Statement

**Academic Reproducibility & Open Science:**
- **Source Code:** The full CuraCascade bioinformatics pipeline and web dashboard are open-source (MIT License).
- **Containerization:** The platform is distributed as a Docker image for isolated, reproducible execution across clinical and academic environments.
- **Data Assets:** The `curacascade_full.db` SQLite/DuckDB relational warehouse and pre-computed `audit_findings.csv` exports are deposited in Zenodo (DOI: 10.5281/zenodo.XXXXXXX) ensuring long-term forensic provenance.

---

## 📚 Methodological References

1. **Richards, S., et al. (2015).** Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine*, 17(5), 405-424.
2. **Rehm, H. L., et al. (2015).** ClinGen—the Clinical Genome Resource. *New England Journal of Medicine*, 372(23), 2235-2242.
3. **Landrum, M. J., et al. (2020).** ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Research*, 48(D1), D835-D844.
4. **Martin, A. R., et al. (2019).** PanelApp crowdsources expert knowledge to establish consensus diagnostic gene panels. *Nature Genetics*, 51(11), 1560-1565.
